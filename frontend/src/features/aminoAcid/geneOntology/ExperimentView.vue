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
          <div :class="tiles.one.current"
               style="transition: all .1s linear;">
            <div class="q-ma-sm">
              <q-btn size="xs" flat color="positive" style="width: 100%" class="q-mb-xs" label="EXPAND"
                     @click="expandTile('one');"/>
              <AminoAcidTable :on-amino-acid-open="setActiveAminoAcid" :rows="experiment?.aminoAcids"/>
            </div>
          </div>
          <div :class="tiles.two.current"
               style="transition: all .1s linear;">
            <div class="q-ma-sm">
              <q-btn size="xs" flat color="positive" style="width: 100%" class="q-mb-xs" label="EXPAND"
                     @click="expandTile('two');"/>
              <GeneOntologyTree :obo-graph="activeAminoAcid?.go" :key="activeAminoAcid?.name"/>
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
import {defineComponent} from 'vue'
import {QSpinnerOrbit, QVueGlobals} from 'quasar';
import ExperimentHeader from "src/components/ExperimentHeader.vue";
import {Experiment} from "src/features/aminoAcid/types";
import AminoAcidInferenceForm from "src/features/aminoAcid/AminoAcidInferenceForm.vue";
import GeneOntologyTree from "src/features/aminoAcid/geneOntology/GeneOntologyTree.vue";
import {AminoAcid} from "src/features/aminoAcid/geneOntology/types";
import useGeneOntologyStore from "src/features/aminoAcid/geneOntology/storage";
import AminoAcidTable from "src/features/aminoAcid/AminoAcidTable.vue";


export default defineComponent({
  name: 'LocalisationExperimentView',
  data() {
    const store = useGeneOntologyStore();

    return {
      tiles: {
        one: {
          current: 'col-md-3',
          hover: 'col-md-6',
          leave: 'col-md-3',
          otherHover: 'col-md-1'
        },
        two: {
          current: 'col-md-9',
          hover: 'col-md-11',
          leave: 'col-md-9',
          otherHover: 'col-md-6'
        }
      } as { [index: string]: { current: String, hover: String, leave: String, otherHover: String } },
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
    expandTile(index: string) {
      if (this.tiles[index].current === this.tiles[index].hover) {
        for (const key in this.tiles) {
          this.tiles[key].current = this.tiles[key].leave;
        }
        return;
      }

      for (const key in this.tiles) {
        if (key == index) {
          this.tiles[key].current = this.tiles[key].hover;
        } else {
          this.tiles[key].current = this.tiles[key].otherHover;
        }
      }
    },
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
    GeneOntologyTree,
    ExperimentHeader,
    AminoAcidInferenceForm
  }
})
</script>
  