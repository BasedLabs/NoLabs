<script lang="ts">
import {defineComponent} from "vue";
import PdbViewer from '../../components/PdbViewer.vue';
import useSmallMoleculesDesignStore from "./storage";
import {QSpinnerOrbit} from "quasar";
import LogsModal from "./components/LogsModal.vue";
import ExperimentControlButtons from "./components/ExperimentControlButtons.vue";
import ExperimentHeader from "../../components/ExperimentHeader.vue";

interface Data {
  smilesData: Smiles[];
  experiment: Experiment | null;
  interval: any;
  logsModalVisible: boolean;
}

export default defineComponent({
  name: "ExperimentView",
  components: {ExperimentHeader, ExperimentControlButtons, LogsModal, PdbViewer},
  store: useSmallMoleculesDesignStore(),
  experimentId: null as unknown as string,
  smilesColumns: [
    {
      name: 'smiles',
      label: 'Formula',
      align: 'left',
      sortable: false
    },
    {
      name: 'drugLikeness',
      label: 'Drug Likeness',
      align: 'left',
      sortable: true
    },
    {
      name: 'score',
      label: 'Score',
      align: 'left',
      sortable: true
    },
    {
      name: 'createdAt',
      label: 'Generated At',
      align: 'left',
      sortable: true
    }
  ],
  data(): Data {
    return {
      smilesData: [] as Smiles[],
      experiment: null as Experiment | null,
      interval: null,
      logsModalVisible: false
    }
  },
  computed: {
    readonly() {
      return this.experiment!.running || this.experiment!.learningCompleted;
    },
    smilesDataLength() {
      return this.smilesData.length;
    }
  },
  methods: {
    async submit() {
      await this.loader('Saving...', async () => {
        await this.$options.store.saveProperties(this.$options.experimentId, this.experiment?.properties);
      })
    },
    toggleLogs() {
      this.logsModalVisible = !this.logsModalVisible;
    },
    async loader(title: string, operation: any) {
      try {
        this.$q.loading.show({
          spinner: QSpinnerOrbit,
          message: title
        });
        await operation();
      } finally {
        this.$q.loading.hide();
      }
    },
    async startExperiment(id: string) {
      await this.loader('Starting experiment', async () => {
        await this.$options.store.startExperiment(this.$options.experimentId, id);
        this.experiment!.running = true;
      });
    },
    async stopExperiment(id: string) {
      this.$q.dialog({
        title: 'Confirm',
        message: 'Would you like to stop the experiment? The progress will be lost.',
        cancel: true,
        persistent: true
      }).onOk(async () => {
        await this.loader('Stopping experiment', async () => {
          await this.$options.store.stopJob(this.$options.experimentId, id);
          this.experiment!.running = false;
        });
      });
    },
    async generateMoreMolecules() {
      await this.loader('Starting job', async () => {

      });
    },
    async resumeExperiment() {

    },
    async onExperimentNameChange(newName: string) {
      await this.loader('Renaming', async () => {
        await this.$options.store.changeExperimentName(this.$options.experimentId, newName);
      });
    }
  },
  async mounted() {
    await this.loader('Loading experiment', async () => {
      this.$options.experimentId = this.$route.params.experimentId as string;
      this.experiment = await this.$options.store.getExperiment(this.$options.experimentId);
      console.log(this.experiment);
      this.interval = setInterval(async () => {
        this.smilesData = await this.$options.store.smilesData(this.$options.experimentId);
      }, 1000);
    })
  },
  unmounted() {
    if (this.interval) {
      clearInterval(this.interval);
    }
  }
})
</script>

<template>
  <ExperimentHeader v-if="experiment" :experiment-name="experiment?.name"
                    :on-experiment-name-change-submit="onExperimentNameChange">
  </ExperimentHeader>
  <q-separator></q-separator>
  <div class="row justify-end q-pa-md" v-if="experiment">
    <div class="q-mx-md" v-if="experiment.running">
      <q-spinner-orbit
        color="white"
        size="xl"
      />
      <q-tooltip :offset="[0, 8]">Loading</q-tooltip>
    </div>
    <q-btn outline align="between" class="q-mx-md" size="lg" color="info" @click="toggleLogs" icon="timeline">Open logs</q-btn>
    <ExperimentControlButtons :resume-experiment="resumeExperiment" :stop-experiment="stopExperiment"
                              :start-experiment="startExperiment" :experiment="experiment"/>
  </div>
  <q-separator></q-separator>
  <div class="row q-col-gutter-md q-ma-sm">
    <div class="col-xs-12 col-md-4">
      <q-table
        title="Generated smiles"
        :rows="smilesData"
        :columns="$options.smilesColumns"
        row-key="smiles"
        :key="smilesDataLength"
      />
    </div>
    <div class="col-xs-12 col-md-5">
      <p class="text-subtitle1">PDB Viewer</p>
      <PdbViewer v-if="experiment?.properties.pdbFile" :pdb-file="experiment!.properties.pdbFile"
                 :key="experiment!.properties.pdbFile?.name"/>
      <div v-if="!experiment?.properties.pdbFile" class="justify-center">
        <q-icon name="warning" color="warning" size="3rem" />
        PDB file is not uploaded
      </div>
    </div>
    <div class="col-xs-12 col-md-3" v-if="experiment" style="max-width: 95%">

      <q-form v-if="experiment!.properties" class="q-gutter-md" @submit="submit">
        <q-scroll-area style="max-width: 100%; height: 70vh" visible>
          <q-file v-if="experiment!.properties" filled bottom-slots accept=".pdb"
                  v-model="experiment!.properties.pdbFile"
                  :readonly="readonly"
                  label=".pdb file" counter>
            <template v-slot:prepend>
              <q-icon name="cloud_upload" @click.stop.prevent/>
            </template>
            <template v-slot:append>
              <q-icon name="close" class="cursor-pointer"/>
            </template>
            <template v-slot:hint>
              Upload .pdb file
            </template>
          </q-file>
          <p>Search space</p>
          <q-input filled v-model="experiment!.properties.centerX" label="Center X" type="number" step="0.1"
                   :readonly="readonly"
                   :rules="[val => val && val > 0 || 'Please type something']">
            <q-tooltip class="text-body1" :offset="[10, 10]" max-width="500px">
              X coordinate of the center of space box
            </q-tooltip>
          </q-input>
          <q-input filled v-model="experiment!.properties.centerY" label="Center Y" type="number" step="0.1"
                   :readonly="readonly"
                   :rules="[val => val && val > 0 || 'Please type something']">
            <q-tooltip class="text-body1" :offset="[10, 10]" max-width="500px">
              Y coordinate of the center of space box
            </q-tooltip>
          </q-input>
          <q-input filled v-model="experiment!.properties.centerZ" label="Center Z" type="number" step="0.1"
                   :readonly="readonly"
                   :rules="[val => val && val > 0 || 'Please type something']">
            <q-tooltip class="text-body1" :offset="[10, 10]" max-width="500px">
              Z coordinate of the center of space box
            </q-tooltip>
          </q-input>
          <q-input filled v-model="experiment!.properties.sizeX" label="Size X" type="number" step="0.1"
                   :readonly="readonly"
                   :rules="[val => val && val > 0 || 'Please type something']">
            <q-tooltip class="text-body1" :offset="[10, 10]" max-width="500px">
              X coordinate of the corner of space box
            </q-tooltip>
          </q-input>
          <q-input filled v-model="experiment!.properties.sizeY" label="Size Y" type="number" step="0.1"
                   :readonly="readonly"
                   :rules="[val => val && val > 0 || 'Please type something']">
            <q-tooltip class="text-body1" :offset="[10, 10]" max-width="500px">
              Y coordinate of the corner of space box
            </q-tooltip>
          </q-input>
          <q-input filled v-model="experiment!.properties.sizeZ" label="Size Z" type="number" step="0.1"
                   :readonly="readonly"
                   :rules="[val => val && val > 0 || 'Please type something']">
            <q-tooltip class="text-body1" :offset="[10, 10]" max-width="500px">
              Z coordinate of the corner of space box
            </q-tooltip>
          </q-input>
          <p>Other parameters</p>
          <p>1 epoch training will take 20-30 min approx</p>
          <q-input filled v-model="experiment!.properties.batchSize" label="Epoch batch size" type="number" step="1"
                   :readonly="readonly"
                   :rules="[val => val && val > 0 || 'Please type something']">
          </q-input>
          <q-input filled v-model="experiment!.properties.minscore" label="Minimum docking acceptance score"
                   type="number" step="0.01"
                   :readonly="readonly"
                   :rules="[val => val && val > 0 || 'Please type something']">
          </q-input>
          <q-input filled v-model="experiment!.properties.epochs" label="Epochs" type="number" step="1"
                   :readonly="readonly"
                   :rules="[val => val && val > 0 || 'Please type something']">
          </q-input>
        </q-scroll-area>
        <q-btn outline style="width: 95%;" label="Save parameters" size="md" type="submit" color="info" :disabled="readonly"/>
      </q-form>
    </div>
  </div>
  <LogsModal :experiment-id="$options.experimentId" v-model:visible="logsModalVisible"/>
</template>
