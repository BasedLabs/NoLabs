<script lang="ts">
import {defineComponent, toRaw} from "vue";
import PdbViewer from '../../components/PdbViewer.vue';
import useSmallMoleculesDesignStore from "./storage";
import {exportFile, Notify, QSpinnerOrbit} from "quasar";
import LogsModal from "./components/LogsModal.vue";
import ExperimentControlButtons from "./components/ExperimentControlButtons.vue";
import ExperimentHeader from "../../components/ExperimentHeader.vue";

interface ShapeCoordinates {
  x: number;
  y: number;
  z: number;
  sizeX: number;
  sizeY: number;
  sizeZ: number;
}

interface Data {
  smilesData: Smiles[];
  experiment: Experiment | null;
  interval: any;
  logsModalVisible: boolean;
  pdbShape: any,
  pdbStage: any,
  shapeCoordinates: ShapeCoordinates | null
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
      sortable: false,
      required: true,
      field: 'smiles'
    },
    {
      name: 'drugLikeness',
      label: 'Drug Likeness',
      align: 'left',
      sortable: true,
      required: true,
      field: 'drugLikeness'
    },
    {
      name: 'score',
      label: 'Binding Score',
      align: 'left',
      sortable: true,
      required: true,
      field: 'score'
    },
    {
      name: 'stage',
      label: 'Stage',
      align: 'left',
      sortable: true,
      required: true,
      field: 'stage'
    }
  ],
  data(): Data {
    return {
      smilesData: [] as Smiles[],
      experiment: null as Experiment | null,
      interval: null,
      logsModalVisible: false,
      pdbShape: null,
      pdbStage: null,
      shapeCoordinates: null
    }
  },
  computed: {
    readonly() {
      return this.experiment!.running;
    },
  },
  watch: {
    shapeCoordinates(newCoordinates: ShapeCoordinates) {
      if (this.pdbStage) {
        setTimeout(() => {
          this.renderSearchBox(
            this.pdbStage,
            newCoordinates.x,
            newCoordinates.y,
            newCoordinates.z,
            newCoordinates.sizeX,
            newCoordinates.sizeY,
            newCoordinates.sizeZ
          );
        }, 100);
      }
    },
    experiment: {
      handler(newExperiment: Experiment) {
        if (!this.shapeCoordinates ||
          newExperiment.properties.centerX != this.shapeCoordinates.x ||
          newExperiment.properties.centerY != this.shapeCoordinates.y ||
          newExperiment.properties.centerZ != this.shapeCoordinates.z ||
          newExperiment.properties.sizeX != this.shapeCoordinates.sizeX ||
          newExperiment.properties.sizeY != this.shapeCoordinates.sizeY ||
          newExperiment.properties.sizeZ != this.shapeCoordinates.sizeZ) {
          this.shapeCoordinates = {
            x: newExperiment.properties.centerX,
            y: newExperiment.properties.centerY,
            z: newExperiment.properties.centerZ,
            sizeX: newExperiment.properties.sizeX,
            sizeY: newExperiment.properties.sizeY,
            sizeZ: newExperiment.properties.sizeZ
          }
        }
      }, deep: true
    }
  },
  methods: {
    copyContent(s: string) {
      navigator.clipboard.writeText(s);
      Notify.create({
        type: "positive",
        closeBtn: 'Close',
        message: `Copied ${s}`
      });
    },
    onRender(stage: any) {
      setTimeout(() => {
        this.renderSearchBox(
          stage,
          this.experiment!.properties.centerX,
          this.experiment!.properties.centerY,
          this.experiment!.properties.centerZ,
          this.experiment!.properties.sizeX,
          this.experiment!.properties.sizeY,
          this.experiment!.properties.sizeZ);

        this.pdbStage = stage;
      }, 1000);
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
    renderSearchBox(stage: any, x: number, y: number, z: number, sizeX: number, sizeY: number, sizeZ: number) {
      if (this.pdbShape) {
        try {
          toRaw(this.pdbShape).removeAllRepresentations();
          stage.removeComponent(toRaw(this.pdbShape));
        } catch {
        }
      }

      const shape = new NGL.Shape("shape");

      var boxBuffer = new NGL.BoxBuffer({
        position: new Float32Array([x, y, z]),
        color: new Float32Array([1, 0, 1]),
        size: new Float32Array([sizeX]),
        heightAxis: new Float32Array([0, sizeY, 0]),
        depthAxis: new Float32Array([0, 0, sizeZ])
      });
      shape.addBuffer(boxBuffer)
      const shapeComp = toRaw(stage).addComponentFromObject(shape)
      shapeComp.addRepresentation("buffer", {opacity: 0.4});
      this.pdbShape = shapeComp;
    },
    onAtomClick(stage: any, coordinates: { x: number, y: number, z: number }) {
      this.pdbStage = stage;
      this.experiment!.properties.centerX = coordinates.x;
      this.experiment!.properties.centerY = coordinates.y;
      this.experiment!.properties.centerZ = coordinates.z;
    },
    onPdbInputClear() {
      this.experiment!.properties.pdbFile = null;
    },
    exportLigandsToCsv() {
      const wrapCsvValue = (val: any, formatFn: any, row: any) => {
        let formatted = formatFn !== void 0
          ? formatFn(val, row)
          : val

        formatted = formatted === void 0 || formatted === null
          ? ''
          : String(formatted)

        formatted = formatted.split('"').join('""')

        return `"${formatted}"`
      }

      const content = [this.$options.smilesColumns.map(col => wrapCsvValue(col.label))].concat(
        this.smilesData.map(row => this.$options.smilesColumns.map(col => wrapCsvValue(
          typeof col.field === 'function'
            ? col.field(row)
            : row[col.field === void 0 ? col.name : col.field],
          col.format,
          row
        )).join(','))
      ).join('\r\n')

      const status = exportFile(
        'table-export.csv',
        content,
        'text/csv'
      )

      if (status !== true) {
        Notify.create({
          message: 'Browser denied file download...',
          color: 'negative',
          icon: 'warning'
        })
      }
    },
    validateProperties(): boolean {
      if (!this.experiment?.properties.pdbFile) {
        Notify.create({
          type: "negative",
          closeBtn: 'Close',
          message: 'Pdb file is empty'
        });
        return false;
      }

      if (this.experiment?.properties.centerX === 0 && this.experiment?.properties.centerY === 0 && this.experiment?.properties.centerZ === 0) {
        Notify.create({
          type: "negative",
          closeBtn: 'Close',
          message: 'Search space center has a X,Y,Z = 0'
        });
        return false;
      }

      if (this.experiment?.properties.sizeX === 0 || this.experiment?.properties.sizeY === 0 || this.experiment?.properties.sizeZ === 0) {
        Notify.create({
          type: "negative",
          closeBtn: 'Close',
          message: 'Size of search box is 0'
        });
        return false;
      }

      return true;
    },
    async startExperiment(id: string) {
      if (!this.validateProperties()) {
        return;
      }

      this.$q.dialog({
        title: 'Confirm',
        message: "Are you sure that you want to start AI learning with following properties? This will overwrite all existing data",
        cancel: true,
        persistent: true,

      }).onOk(async () => {
        await this.loader('Saving properties', async () => {
          await this.$options.store.saveProperties(this.$options.experimentId, this.experiment?.properties);
        });

        await this.loader('Starting experiment', async () => {
          await this.$options.store.startExperiment(this.$options.experimentId, id);
          this.experiment!.running = true;
        });
      });
    },
    async stopExperiment() {
      this.$q.dialog({
        title: 'Confirm',
        message: 'Would you like to stop the experiment? You will be able to run sampling, but learning will be frozen on the last stage',
        cancel: true,
        persistent: true
      }).onOk(async () => {
        await this.loader('Stopping experiment', async () => {
          await this.$options.store.stopExperiment(this.$options.experimentId);
          this.experiment!.running = false;
        });
      });
    },
    async startSampling() {
      this.$q.dialog({
        title: 'Sampling',
        message: 'Enter number of molecules to generate',
        prompt: {
          model: '1',
          type: 'number' // optional
        },
        cancel: true,
        persistent: true
      }).onOk(async data => {
        await this.loader('Starting sampling', async () => {
          await this.$options.store.startSampling(this.$options.experimentId, data);
          this.experiment!.running = true;
        });
      });
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
      this.smilesData = await this.$options.store.smilesData(this.$options.experimentId);
      this.interval = setInterval(async () => {
        this.smilesData = await this.$options.store.smilesData(this.$options.experimentId);
        const status = await this.$options.store.status(this.$options.experimentId);
        this.experiment!.running = status.running;
        this.experiment!.samplingAllowed = status.samplingAllowed;
      }, 5000);
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
    <q-banner class="bg-primary text-white">
      <q-badge color="positive">
        HINT
      </q-badge>
      Click on the atom to fill-up the center of the search box. Enter SizeXYZ to set size of the search box.
    </q-banner>
    <div class="q-mx-md" v-if="experiment.running">
      <q-spinner-orbit
        color="white"
        size="xl"
      />
      <q-tooltip :offset="[0, 8]">Loading</q-tooltip>
    </div>
    <q-btn outline align="between" class="q-mx-md" size="md" color="info" @click="toggleLogs" icon="timeline">Open
      logs
    </q-btn>
    <ExperimentControlButtons :start-sampling="startSampling" :stop-experiment="stopExperiment"
                              :start-experiment="startExperiment" :experiment="experiment"/>
  </div>
  <q-separator></q-separator>
  <div class="row q-col-gutter-md q-ma-sm">
    <div class="col-xs-12 col-md-5">
      <q-table
        title="Generated smiles"
        :rows="smilesData"
        :columns="$options.smilesColumns"
        separator="cell"
      >
        <template v-slot:body="props">
          <q-tr :props="props">
            <q-td key="smiles" :props="props" class="row justify-center">
              <q-btn outline color="info" icon="content_copy" size="xs" @click="copyContent(props.row.smiles)"/>
              <q-input class="q-pl-xs" v-model="props.row.smiles" readonly dense/>
            </q-td>
            <q-td key="drugLikeness" :props="props">
              {{ props.row.drugLikeness.toFixed(2) }}
            </q-td>
            <q-td key="score" :props="props">
              {{ props.row.score.toFixed(2) }}
            </q-td>
            <q-td key="stage" :props="props">
              {{ props.row.stage }}
            </q-td>
          </q-tr>
        </template>
        <template v-slot:top-right>
          <q-btn
            color="positive"
            outline
            icon-right="archive"
            label="Export to csv"
            no-caps
            @click="exportLigandsToCsv"
          />
        </template>
      </q-table>
    </div>
    <div class="col-xs-12 col-md-5">
      <p class="text-subtitle1">PDB Viewer</p>
      <PdbViewer v-if="experiment?.properties.pdbFile" :pdb-file="experiment!.properties.pdbFile"
                 :on-atom-click="onAtomClick"
                 :on-render="onRender"
                 :key="experiment!.properties.pdbFile?.name"/>
      <div v-if="!experiment?.properties.pdbFile" class="justify-center">
        <q-icon name="warning" color="warning" size="3rem"/>
        PDB file is not uploaded
      </div>
    </div>
    <div class="col-xs-12 col-md-2" v-if="experiment" style="max-width: 95%">

      <q-form v-if="experiment!.properties" class="q-gutter-md">
        <q-scroll-area style="max-width: 100%; height: 70vh" visible>
          <q-file v-if="experiment!.properties" filled bottom-slots accept=".pdb"
                  v-model="experiment!.properties.pdbFile"
                  :readonly="readonly"
                  clearable
                  @clear="onPdbInputClear"
                  label=".pdb file" counter>
            <template v-slot:prepend>
              <q-icon name="cloud_upload" @click.stop.prevent/>
            </template>
            <template v-slot:hint>
              Upload .pdb file
            </template>
          </q-file>
          <p>Search space</p>
          <q-input filled v-model="experiment!.properties.centerX" label="Center X" type="number" step="0.01"
                   :readonly="readonly"
                   :rules="[val => val !== undefined && val !== null || 'Invalid input']">
            <q-tooltip class="text-body1" :offset="[10, 10]" max-width="500px">
              X coordinate of the center of space box
            </q-tooltip>
          </q-input>
          <q-input filled v-model="experiment!.properties.centerY" label="Center Y" type="number" step="0.01"
                   :readonly="readonly"
                   :rules="[val => val !== undefined && val !== null || 'Invalid input']">
            <q-tooltip class="text-body1" :offset="[10, 10]" max-width="500px">
              Y coordinate of the center of space box
            </q-tooltip>
          </q-input>
          <q-input filled v-model="experiment!.properties.centerZ" label="Center Z" type="number" step="0.01"
                   :readonly="readonly"
                   :rules="[val => val !== undefined && val !== null || 'Invalid input']">
            <q-tooltip class="text-body1" :offset="[10, 10]" max-width="500px">
              Z coordinate of the center of space box
            </q-tooltip>
          </q-input>
          <q-input filled v-model="experiment!.properties.sizeX" label="Size X" type="number" step="0.01"
                   :readonly="readonly"
                   :rules="[val => val !== undefined && val !== null && val > 0 && val <= 30 || 'Invalid input']">
            <q-tooltip class="text-body1" :offset="[10, 10]" max-width="500px">
              X coordinate of the corner of space box
            </q-tooltip>
          </q-input>
          <q-input filled v-model="experiment!.properties.sizeY" label="Size Y" type="number" step="0.01"
                   :readonly="readonly"
                   :rules="[val => val !== undefined && val !== null && val > 0 && val <= 30 || 'Invalid input']">
            <q-tooltip class="text-body1" :offset="[10, 10]" max-width="500px">
              Y coordinate of the corner of space box
            </q-tooltip>
          </q-input>
          <q-input filled v-model="experiment!.properties.sizeZ" label="Size Z" type="number" step="0.01"
                   :readonly="readonly"
                   :rules="[val => val !== undefined && val !== null && val > 0 && val <= 30 || 'Invalid input']">
            <q-tooltip class="text-body1" :offset="[10, 10]" max-width="500px">
              Z coordinate of the corner of space box
            </q-tooltip>
          </q-input>
          <p>Other parameters</p>
          <p>1 epoch training will take 20-30 min approx</p>
          <q-input filled v-model="experiment!.properties.batchSize" label="Epoch batch size" type="number" step="1"
                   :readonly="readonly"
                   :rules="[val => val !== undefined && val !== null && val > 0 || 'Invalid input']">
          </q-input>
          <q-input filled v-model="experiment!.properties.minscore" label="Minimum docking acceptance score"
                   type="number" step="0.01"
                   :readonly="readonly"
                   :rules="[val => val !== undefined && val !== null || 'Invalid input']">
          </q-input>
          <q-input filled v-model="experiment!.properties.epochs" label="Epochs" type="number" step="1"
                   :readonly="readonly"
                   :rules="[val => val !== undefined && val !== null && val > 0 || 'Invalid input']">
          </q-input>
        </q-scroll-area>
      </q-form>
    </div>
  </div>
  <LogsModal :experiment-id="$options.experimentId" v-model:visible="logsModalVisible"/>
</template>
