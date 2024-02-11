<template>
  <div v-if="experimentLoaded">
    <q-separator></q-separator>
    <q-layout container style="height: 100vh">
      <ExperimentHeader :experiment-name="experiment?.name" :on-experiment-name-change-submit="onExperimentNameChange">
        <q-btn color="positive" size="md" outline label="Simulation parameters"
               @click="showInferenceForm = !showInferenceForm"/>
      </ExperimentHeader>
      <q-page-container>
        <div class="row">
          <div :class="tiles.one.current"
               style="transition: all .1s linear;">
            <div class="q-ma-sm">
              <q-btn size="xs" flat color="positive" style="width: 100%" class="q-mb-xs" label="EXPAND"
                     @click="expandTile('one');"/>
              <PdbViewer v-if="experimentHasInputData" :pdb-file="experiment?.properties.pdbFile"
                         :key="experiment?.pdbContent?.name"/>
            </div>
          </div>
          <div :class="tiles.two.current"
               style="transition: all .1s linear;">
            <div class="q-pl-sm q-ma-sm">
              <q-btn size="xs" flat color="positive" style="width: 100%" class="q-mb-xs" label="EXPAND"
                     @click="expandTile('two');"/>
              <TimelineView :timeline="experiment?.timeline"/>
            </div>
          </div>
          <div :class="tiles.three.current" style="transition: all .1s linear;">
            <div class="q-mt-sm q-mb-sm q-mr-sm">
              <q-btn flat size="xs" color="positive" style="width: 100%" class="q-mb-xs" label="EXPAND"
                     @click="expandTile('three');"/>
              <h6 v-if="simulationCriticalError" class="q-ma-md">
                <q-icon name="warning" color="warning" size="4rem"/>
                Simulations error. Fix errors (check them in table) and restart simulation.
              </h6>
              <PdbViewer v-if="!simulationCriticalError && experimentHasGeneratedData"
                         :key="experiment?.pdbContent?.name"
                         :pdb-file="experiment?.pdbContent" :simulation="true"
                         file-name-prefix="Simulation"/>
            </div>
          </div>
        </div>
      </q-page-container>
    </q-layout>
    <q-dialog v-model="showInferenceForm" position="right" :maximized="true">
      <q-card>
        <q-card-section class="row items-center q-pb-none">
          <div class="text-h6">Simulation parameters</div>
          <q-space/>
          <q-btn icon="close" flat round dense v-close-popup/>
        </q-card-section>
        <q-card-section>
          <InferenceFormView :on-submit="onSubmit" :properties="experiment?.properties"/>
        </q-card-section>
      </q-card>
    </q-dialog>
  </div>
</template>

<script lang="ts">
import {defineComponent} from 'vue'
import {QVueGlobals, QSpinnerOrbit} from 'quasar';
import PdbViewer from "src/components/PdbViewer.vue";
import InferenceFormView from "src/features/conformations/InferenceFormView.vue";
import useConformationsStore from "src/features/conformations/storage";
import {Experiment, ExperimentProperties} from "src/features/conformations/types";
import TimelineView from "src/components/Timeline.vue";
import ExperimentHeader from "src/components/ExperimentHeader.vue";


export default defineComponent({
  name: 'ConformationExperimentView',
  data() {
    const store = useConformationsStore();

    return {
      experiment: null as Experiment,
      showInferenceForm: false,
      store,
      tiles: {
        one: {
          current: 'col-4',
          hover: 'col-6',
          leave: 'col-4',
          otherHover: 'col-3'
        },
        two: {
          current: 'col-4',
          hover: 'col-6',
          leave: 'col-4',
          otherHover: 'col-3'
        },
        three: {
          current: 'col-4',
          hover: 'col-6',
          leave: 'col-4',
          otherHover: 'col-3'
        }
      } as { [index: string]: { current: String, hover: String, leave: String, otherHover: String } }
    }
  },
  computed: {
    experimentLoaded(): boolean {
      return this.experiment !== null;
    },
    experimentHasInputData(): boolean {
      return this.experimentLoaded && this.experiment!.properties.pdbFile != null;
    },
    experimentHasGeneratedData(): boolean {
      return this.experimentLoaded && this.experiment!.pdbContent != null && this.experiment!.pdbContent!.size > 4;
    },
    simulationCriticalError(): boolean {
      return this.experimentLoaded && this.experiment?.timeline.length > 0 && (this.experiment?.pdbContent == null || this.experiment?.pdbContent.size <= 4);
    }
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
    async onSubmit(properties: ExperimentProperties) {
      this.$q.loading.show({
        spinner: QSpinnerOrbit,
        message: 'Running computations'
      });

      const response = await this.store.inference({
        experimentId: this.experiment!.id,
        experimentName: this.experiment!.name,
        pdbFile: properties.pdbFile!,
        totalFrames: properties.totalFrames!,
        temperatureK: properties.temperatureK!,
        takeFrameEvery: properties.takeFrameEvery!,
        stepSize: properties.stepSize!,
        replaceNonStandardResidues: properties.replaceNonStandardResidues!,
        addMissingAtoms: properties.addMissingAtoms!,
        addMissingHydrogens: properties.addMissingHydrogens!,
        frictionCoeff: properties.frictionCoeff!,
        ignoreMissingAtoms: properties.ignoreMissingAtoms!,
        integrator: properties.integrator!,
      });

      if (response.experiment !== null) {
        this.experiment = response.experiment;
      }

      this.showInferenceForm = false;

      this.$q.loading.hide();
    },
    async onExperimentNameChange(newExperimentName: string) {
      await this.store.changeExperimentName(this.experiment?.id as string, newExperimentName);
      this.experiment!.name = newExperimentName;
    }
  },
  async mounted() {
    const experimentId = this.$route.params.experimentId as string;

    this.$q.loading.show({
      spinner: QSpinnerOrbit,
      message: `Experiment ${experimentId}`
    });

    const response = await this.store.getExperiment(experimentId);

    if (response.experiment !== null) {
      this.experiment = response.experiment;
    }

    this.$q.loading.hide();

    if (!this.experimentHasInputData) {
      this.showInferenceForm = true;
    }
  },
  components: {
    ExperimentHeader,
    TimelineView,
    PdbViewer,
    InferenceFormView
  }
})
</script>
  