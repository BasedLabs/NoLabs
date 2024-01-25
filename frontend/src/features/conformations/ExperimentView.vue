<template>
  <div v-if="experimentLoaded">
    <q-separator></q-separator>
    <q-layout container style="height: 100vh">
      <q-header :class="$q.dark.isActive ? 'bg-secondary' : 'bg-black'">
        <q-toolbar>
          <q-toolbar-title>{{ experiment.name }}
            <q-btn round
                   @click="changeExperimentName" color="positive" size="sm" flat icon="edit"/>
          </q-toolbar-title>
          <q-btn color="positive" size="md" outline label="Binder parameters"
                 @click="showInferenceForm = !showInferenceForm"/>
        </q-toolbar>
      </q-header>
      <q-page-container>
        <div class="row" v-if="experimentHasGeneratedData">
          <div class="col-6">
            <div class="q-ma-sm">
              <PdbViewer :pdb-file="experiment?.properties.pdbFile"/>
            </div>
          </div>
          <div class="col-6">
            <div class="q-mt-sm q-mb-sm q-mr-sm">
              <PdbViewer :pdb-file="experiment?.pdbContent"/>
            </div>
          </div>
        </div>
      </q-page-container>
    </q-layout>
    <q-dialog v-model="showInferenceForm" position="right" :maximized="true">
      <q-card>
        <q-card-section class="row items-center q-pb-none">
          <div class="text-h6">Protein binder design parameters</div>
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
import {QVueGlobals, useQuasar, QSpinnerOrbit} from 'quasar';
import PdbViewer from "src/components/PdbViewer.vue";
import InferenceFormView from "src/features/conformations/InferenceFormView.vue";
import useConformationsStorage from "src/features/conformations/storage";
import {Experiment, ExperimentProperties} from "src/features/conformations/types";
import type {IntegratorsRequest} from "src/api/client";


export default defineComponent({
  name: 'ProteinDesignExperimentView',
  data() {
    const store = useConformationsStorage();

    return {
      experiment: null as Experiment,
      showInferenceForm: false,
      store,
      quasar: null as unknown as QVueGlobals,
    }
  },
  computed: {
    experimentLoaded(): boolean {
      return this.experiment !== null;
    },
    experimentHasGeneratedData(): boolean {
      return this.experimentLoaded && this.experiment!.pdbContent != null;
    }
  },
  methods: {
    async onSubmit(properties: ExperimentProperties) {
      this.quasar.loading.show({
        spinner: QSpinnerOrbit,
        message: 'Running AI models. This can take a couple of minutes'
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

      this.quasar.loading.hide();
    },
    changeExperimentName() {
      this.quasar.dialog({
        color: 'positive',
        title: 'Prompt',
        message: 'Enter new experiment name',
        prompt: {
          model: this.experiment!.name,
          required: true,
          type: 'text' // optional
        },
        cancel: true,
        persistent: true
      }).onOk(async data => {
        if (!data)
          return;
        this.quasar.loading.show({
          spinner: QSpinnerOrbit,
          message: 'Changing experiment name'
        });
        await this.store.changeExperimentName(this.experiment?.id as string, data);
        this.experiment!.name = data;
        this.quasar.loading.hide();
      });
    }
  },
  async mounted() {
    const experimentId = this.$route.params.experimentId as string;

    this.quasar = useQuasar();

    this.quasar.loading.show({
      spinner: QSpinnerOrbit,
      message: `Experiment ${experimentId}`
    });

    const response = await this.store.getExperiment(experimentId);

    if (response.experiment !== null) {
      this.experiment = response.experiment;
    }

    this.quasar.loading.hide();

    if (!this.experimentHasGeneratedData) {
      this.showInferenceForm = true;
    }
  },
  components: {
    PdbViewer,
    InferenceFormView
  }
})
</script>
  