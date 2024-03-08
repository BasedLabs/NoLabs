<template>
  <div v-if="experimentLoaded">
    <q-separator></q-separator>
    <q-layout container style="height: 100vh">
      <ExperimentHeader :experiment-name="experiment?.name" :on-experiment-name-change-submit="onExperimentNameChange">
        <q-btn color="info" size="md" outline label="Solubility parameters"
               @click="showInferenceForm = !showInferenceForm"/>
      </ExperimentHeader>
      <q-page-container>
        <div class="row" v-if="experimentHasGeneratedData">
          <div class="col-6">
            <div class="q-ma-sm">
              <AminoAcidTable :on-amino-acid-open="setActiveAminoAcid" :rows="experiment?.aminoAcids"/>
            </div>
          </div>
          <div class="col-1">
          </div>
          <div class="col-4">
            <div class="q-pl-sm q-ma-sm">
              <h6>{{ activeAminoAcid.name }} soluble with {{ activeAminoAcid.soluble_probability * 100.0 }}% probability.</h6>
            </div>
          </div>
        </div>
      </q-page-container>
    </q-layout>
    <q-dialog v-model="showInferenceForm" position="right" :maximized="true">
      <q-card>
        <q-card-section class="row items-center q-pb-none">
          <div class="text-h6">Localisation parameters</div>
          <q-space/>
          <q-btn icon="close" flat round dense v-close-popup/>
        </q-card-section>
        <q-card-section>
          <AminoAcidInferenceForm :on-submit="onSubmit" :properties="experiment?.properties"/>
        </q-card-section>
      </q-card>
    </q-dialog>
  </div>
</template>

<script lang="ts">
import {defineComponent, ref} from 'vue'
import {QVueGlobals, QSpinnerOrbit} from 'quasar';
import ExperimentHeader from "src/components/ExperimentHeader.vue";
import {Experiment} from "src/features/aminoAcid/types";
import {AminoAcid} from "src/features/aminoAcid/solubility/types";
import useSolubilityStore from "src/features/aminoAcid/solubility/storage";
import AminoAcidInferenceForm from "src/features/aminoAcid/AminoAcidInferenceForm.vue";
import AminoAcidTable from "src/features/aminoAcid/AminoAcidTable.vue";


export default defineComponent({
  name: 'LocalisationExperimentView',
  data() {
    const store = useSolubilityStore();

    return {
      experiment: null as Experiment<AminoAcid>,
      showInferenceForm: false,
      store,
      activeAminoAcid: null as AminoAcid | null | undefined
    }
  },
  computed: {
    experimentLoaded(): boolean {
      return this.experiment !== null;
    },
    experimentHasGeneratedData(): boolean {
      return this.experimentLoaded && this.experiment!.aminoAcids.length > 0;
    },
    aminoAcidRows() {
      return this.experiment?.aminoAcids;
    },
  },
  methods: {
    setActiveAminoAcid(aminoAcidName: string): void {
      const aminoAcid = this.experiment?.aminoAcids.find(x => x.name === aminoAcidName);
      this.activeAminoAcid = aminoAcid;
    },
    setExperiment(experiment: Experiment<AminoAcid> | null){
      if (experiment !== null) {
        this.experiment = experiment;

        if(experiment.aminoAcids.length > 0){
          this.setActiveAminoAcid(experiment.aminoAcids[0]!.name);
        }
      }
    },
    async onExperimentNameChange(newExperimentName: string) {
      await this.store.changeExperimentName(this.experiment?.id as string, newExperimentName);
      this.experiment!.name = newExperimentName;
    },
    async onSubmit(data: { aminoAcidSequence: string, fastas: Array<File> }) {
      this.$q.loading.show({
        spinner: QSpinnerOrbit,
        message: 'Running AI models. This can take a couple of minutes'
      });

      const response = await this.store.inference({
        experimentId: this.experiment!.id,
        experimentName: this.experiment!.name,
        aminoAcidSequence: data.aminoAcidSequence,
        fastas: data.fastas
      });

      this.setExperiment(response.experiment);

      this.showInferenceForm = false;

      this.$q.loading.hide();
    },
  },
  async mounted() {
    const experimentId = this.$route.params.experimentId as string;

    this.$q.loading.show({
      spinner: QSpinnerOrbit,
      message: `Experiment ${experimentId}`
    });

    const response = await this.store.getExperiment(experimentId);

    this.setExperiment(response.experiment);

    this.$q.loading.hide();

    if (!this.experimentHasGeneratedData) {
      this.showInferenceForm = true;
    }
  },
  components: {
    AminoAcidTable,
    ExperimentHeader,
    AminoAcidInferenceForm
  }
})
</script>
