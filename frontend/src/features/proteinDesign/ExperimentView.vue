<template>
  <div v-if="experimentLoaded">
    <q-separator></q-separator>
    <q-layout container style="height: 100vh">
      <ExperimentHeader :experiment-name="experiment?.name" :on-experiment-name-change-submit="onExperimentNameChange">
        <q-btn color="info" size="md" outline label="Binder parameters"
               @click="showInferenceForm = !showInferenceForm"/>
      </ExperimentHeader>
      <q-page-container>
        <div class="row" v-if="experimentHasGeneratedData">
          <div class="col-5">
            <div class="q-ma-sm">
              <PdbViewer :pdb-file="experiment?.properties.inputPdbFile"/>
            </div>
          </div>
          <div class="col-2">
            <div class="q-pl-sm q-ma-sm">
              <q-table
                  title="Generated pdbs"
                  :rows="generatedPdbsTableRows"
                  :columns="generatedPdbsTableColumns"
                  row-key="name"
                  @row-click="(_, row) => {generatedPdbsReload = !generatedPdbsReload; selectedGeneratedPdbIndex = row.id;}"
                  :pagination="{rowsPerPage: 5}"
              />
            </div>
          </div>
          <div class="col-5">
            <div class="q-mt-sm q-mb-sm q-mr-sm">
              <PdbViewer :pdb-file="experiment!.generatedPdbs[selectedGeneratedPdbIndex]"
                         :key="generatedPdbsReload"/>
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
import {defineComponent, ref} from 'vue'
import useProteinDesignStore from 'src/features/proteinDesign/storage';
import {QVueGlobals, QSpinnerOrbit} from 'quasar';
import InferenceFormView from "src/features/proteinDesign/InferenceFormView.vue";
import {Experiment, ExperimentProperties} from "src/features/proteinDesign/types";
import PdbViewer from "src/components/PdbViewer.vue";
import ExperimentHeader from "src/components/ExperimentHeader.vue";


export default defineComponent({
  name: 'ProteinDesignExperimentView',
  data() {
    const store = useProteinDesignStore();

    return {
      experiment: null as Experiment,
      showInferenceForm: false,
      store,
      selectedGeneratedPdbIndex: 0,
      generatedPdbsReload: false
    }
  },
  computed: {
    generatedPdbsTableColumns() {
      return [
        {
          name: 'id',
          required: true,
          label: '#',
          align: 'left',
          sortable: false,
          field: row => row.id,
        },
        {
          name: 'name',
          required: true,
          label: 'Name',
          align: 'left',
          sortable: false,
          field: row => row.name
        }
      ]
    },
    generatedPdbsTableRows(): { id: number, name: string }[] {
      return this.experiment ? this.experiment.generatedPdbs.map((file, i) => {
        return {
          id: i, name: file.name
        }
      }) : [];
    },
    experimentLoaded(): boolean {
      return this.experiment !== null;
    },
    experimentHasGeneratedData(): boolean {
      return this.experimentLoaded && this.experiment!.generatedPdbs.length > 0;
    }
  },
  methods: {
    async onExperimentNameChange(newExperimentName: string) {
      await this.store.changeExperimentName(this.experiment?.id as string, newExperimentName);
      this.experiment!.name = newExperimentName;
    },
    async onSubmit(properties: ExperimentProperties) {
      this.$q.loading.show({
        spinner: QSpinnerOrbit,
        message: 'Running AI models. This can take a couple of minutes'
      });

      const response = await this.store.inference({
        experimentId: this.experiment?.id,
        experimentName: this.experiment?.name as string,
        pdbFile: properties.inputPdbFile!,
        contig: properties.contig,
        numberOfDesigns: properties.numberOfDesigns,
        timesteps: properties.timesteps,
        hotspots: properties.hotspots
      });

      if (response.experiment !== null) {
        this.experiment = response.experiment;
      }

      this.showInferenceForm = false;

      this.$q.loading.hide();
    },
    changeExperimentName() {
      this.$q.dialog({
        color: 'info',
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
        this.$q.loading.show({
          spinner: QSpinnerOrbit,
          message: 'Changing experiment name'
        });
        await this.store.changeExperimentName(this.experiment?.id as string, data);
        this.experiment!.name = data;
        this.$q.loading.hide();
      });
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

    if (!this.experimentHasGeneratedData) {
      this.showInferenceForm = true;
    }
  },
  components: {
    ExperimentHeader,
    PdbViewer,
    InferenceFormView
  }
})
</script>
  