<template>
  <div v-if="experimentLoaded">
    <q-separator></q-separator>
    <q-layout container style="height: 100vh">
      <ExperimentHeader :experiment-name="experiment?.name" :on-experiment-name-change-submit="onExperimentNameChange">
        <q-btn color="positive" size="md" outline label="Localisation parameters"
               @click="showInferenceForm = !showInferenceForm"/>
      </ExperimentHeader>
      <q-page-container>
        <div class="row" v-if="experimentHasGeneratedData">
          <div class="col-6">
            <div class="q-ma-sm">
              <AminoAcidTable :on-amino-acid-open="setActiveAminoAcid" :rows="experiment?.aminoAcids"/>
            </div>
          </div>
          <div class="col-6">
            <div class="q-pl-sm q-ma-sm">
              <PdbViewer :pdb-file="activeAminoAcid?.pdbFile" :key="activeAminoAcid?.name"/>
            </div>
          </div>
        </div>
      </q-page-container>
    </q-layout>
    <q-dialog v-model="showInferenceForm" position="right" :maximized="true">
      <q-card>
        <q-card-section class="row items-center q-pb-none">
          <div class="text-h6">Folding parameters</div>
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
import {defineComponent} from 'vue'
import {QSpinnerOrbit, QVueGlobals, useQuasar} from 'quasar';
import ExperimentHeader from "src/components/ExperimentHeader.vue";
import {AminoAcid} from "src/features/aminoAcid/folding/types";
import {Experiment} from "src/features/aminoAcid/types";
import AminoAcidInferenceForm from "src/features/aminoAcid/AminoAcidInferenceForm.vue";
import AminoAcidTable from "src/features/aminoAcid/AminoAcidTable.vue";
import useFoldingStore from "src/features/aminoAcid/folding/storage";
import PdbViewer from "src/components/PdbViewer.vue";


export default defineComponent({
  name: 'LocalisationExperimentView',
  data() {
    const store = useFoldingStore();

    return {
      experiment: null as Experiment<AminoAcid>,
      showInferenceForm: false,
      store,
      quasar: null as unknown as QVueGlobals,
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
      this.activeAminoAcid = this.experiment?.aminoAcids.find(x => x.name === aminoAcidName);
    },
    setExperiment(experiment: Experiment<AminoAcid> | null) {
      if (experiment !== null) {
        this.experiment = experiment;

        if (experiment.aminoAcids.length > 0) {
          this.setActiveAminoAcid(experiment.aminoAcids[0]!.name);
        }
      }
    },
    async onExperimentNameChange(newExperimentName: string) {
      await this.store.changeExperimentName(this.experiment?.id as string, newExperimentName);
      this.experiment!.name = newExperimentName;
    },
    async onSubmit(data: { aminoAcidSequence: string, fastas: Array<File> }) {
      this.quasar.loading.show({
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

      this.quasar.loading.hide();
    },
  },
  async mounted() {
    const experimentId = this.$route.params.experimentId as string;

    this.quasar = useQuasar();

    this.quasar.loading.show({
      spinner: QSpinnerOrbit,
      message: `Experiment ${experimentId}`
    });

    const response = await this.store.getExperiment(experimentId);

    this.setExperiment(response.experiment);

    this.quasar.loading.hide();

    if (!this.experimentHasGeneratedData) {
      this.showInferenceForm = true;
    }
  },
  components: {
    PdbViewer,
    AminoAcidTable,
    ExperimentHeader,
    AminoAcidInferenceForm
  }
})
</script>
  