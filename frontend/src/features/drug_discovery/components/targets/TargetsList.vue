<template>
  <q-list bordered separator class="bg-black">
    <q-item-label header class="text-white">Upload targets</q-item-label>
    <q-item>
      <q-btn size="md" class="full-width q-pm-sm" push color="info" @click="uploadTargetDialog = true">
        + upload targets
        <q-tooltip>Add targets to the experiment</q-tooltip>
      </q-btn>
    </q-item>
    <q-expansion-item
        popup
        expand-separator
        :content-inset-level="1"
        v-for="target in targets"
        :key="target.target_id"
        :label="target.target_name"
        @show="() => selectTarget(target)"
    >
      <q-btn color="negative" @click="deleteTarget(target)" icon="delete" flat>
        Delete target
      </q-btn>
      <TargetDetail v-if="target.data" :experiment-id="experimentId" :original-target="target"></TargetDetail>
    </q-expansion-item>
  </q-list>
  <q-page-sticky position="bottom-right" :offset="[100, 50]">
    <q-btn size="lg" push round color="info" icon="add" @click="uploadTargetDialog = true">
      <q-tooltip>Add targets to the experiment</q-tooltip>
    </q-btn>
  </q-page-sticky>

  <q-dialog v-model="uploadTargetDialog">
    <q-card>
      <q-card-section>
        <div class="text-h6">Upload Target Files</div>
        <q-file v-model="targetFileModel" ref="targetFileInput" accept=".fasta" multiple label="Choose files"/>
      </q-card-section>

      <q-card-actions align="right">
        <q-btn flat label="Close" color="negative" @click="uploadTargetDialog = false"/>
        <q-btn flat label="Upload" color="info" @click="uploadTargetFiles"/>
      </q-card-actions>
    </q-card>
  </q-dialog>
</template>

<script lang="ts">
import TargetDetail from "src/features/drug_discovery/components/targets/TargetDetail.vue";
import {useDrugDiscoveryStore} from "src/features/drug_discovery/storage";
import {QSpinnerOrbit} from "quasar";
import {defineComponent} from "vue";
import {ExtendedTargetMetaData} from "./types";
import {FunctionCallReturnData, RcsbPdbData} from "src/api/client";
import {useBioBuddyStore} from "../../../biobuddy/storage";

export default defineComponent({
  name: "TargetsList",
  components: {
    TargetDetail,
  },
  props: {
    experimentId: {
      type: String,
      required: true,
    },
  },
  data() {
    return {
      uploadTargetDialog: false,
      targetFileModel: [],
      uploadLigandDialog: false,
      uploadingLigandFiles: [],
      selectedTarget: null,
      uploadToAllTargets: false,
      rscbQueryCallBack: null as ((data: FunctionCallReturnData) => void) | null,
      targets: [] as ExtendedTargetMetaData[]
    };
  },
  computed: {
    drugDiscoveryStore() {
      return useDrugDiscoveryStore();
    },
  },
  async mounted() {
    const bioBuddyStore = useBioBuddyStore();
    this.rscbQueryCallBack = async (data: { files: RcsbPdbData[] }) => {
      const files: File[] = [];
      const metaDatasArray: Record<string, string>[] = [];

      for (let rcsbObject of data.files) {
        const file = new File([new Blob([rcsbObject.content!])], "protein.fasta");
        files.push(file);
        const metadataAsString: Record<string, string> = {};
        Object.entries(rcsbObject.metadata).forEach(([key, value]) => {
          metadataAsString[key] = typeof value === 'object' ? JSON.stringify(value) : String(value);
        });
        metaDatasArray.push(metadataAsString);
      }
      await this.handleTargetFileUpload(files, metaDatasArray);
    };
    bioBuddyStore.addQueryRcsbPdbEventHandler(this.rscbQueryCallBack);

    this.targets = await this.drugDiscoveryStore.fetchTargetsForExperiment(this.experimentId) as ExtendedTargetMetaData[];

    },
  methods: {
    async selectTarget(target: ExtendedTargetMetaData) {
      target.data = await this.getTargetData(target);
    },
    uploadTargetFiles() {
      if (this.targetFileModel.length > 0) {
        this.handleTargetFileUpload(this.targetFileModel)
          .then(() => {
            this.targetFileModel = []; // Reset the file input after upload
          })
          .catch(error => {
            console.error("Upload failed", error);
          });
      } else {
        console.log("No files selected");
      }
    },
    async handleTargetFileUpload(files: File[], additionalMetaDataArray?: Record<string, string>[]) {
      const store = useDrugDiscoveryStore();
      this.$q.loading.show({
        spinner: QSpinnerOrbit,
        message: `Uploading target files`
      });

      try {
        for (let index = 0; index < files.length; index++) {
          const file = files[index];
          const metaData = additionalMetaDataArray ? additionalMetaDataArray[index] : undefined;
          const newTargets = await store.uploadTargetToExperiment(this.experimentId, file, metaData);
          newTargets?.forEach(target => this.targets.push(target));
        }
      } catch (error) {
        console.error('Error uploading target files:', error);
      } finally {
        this.uploadTargetDialog = false;
        this.$q.loading.hide();
      }
    },
    async getTargetData(target: ExtendedTargetMetaData) {
      const store = useDrugDiscoveryStore();
      target.loadingTargetData = true;
      try {
        return await store.fetchTargetData(this.experimentId, target.target_id);
      } catch (error) {
        console.error('Error fetching target data:', error);
      } finally {
        target.loadingTargetData = false;
      }
    },
    async deleteTarget(targetToDelete: ExtendedTargetMetaData) {
      const store = useDrugDiscoveryStore();
      this.$q.loading.show({
        spinner: QSpinnerOrbit,
        message: `Deleting target`
      });
      try {
        await store.deleteTargetFromExperiment(this.experimentId, targetToDelete.target_id);
        this.targets = this.targets.filter((target: ExtendedTargetMetaData) => target.target_id !== targetToDelete.target_id);
      } catch (error) {
        console.error('Error deleting target:', error);
      }
      finally{
        this.$q.loading.hide();
      }
    },
  },
  unmounted() {
    const bioBuddyStore = useBioBuddyStore();
    const index = bioBuddyStore.queryRcsbPdbEventHandlers.indexOf(this.rscbQueryCallBack!);
    if (index !== -1) {
      bioBuddyStore.queryRcsbPdbEventHandlers.splice(index, 1);
    }
  }
}
)
</script>
