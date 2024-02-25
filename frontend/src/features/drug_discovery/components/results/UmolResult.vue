<template>
  <div class="umol-result">
    <q-page padding>
      <PdbViewer v-if="predictedPdb && predictedSDF" :pdb-file="predictedPdb" :sdf-file="predictedSDF" />

      <!-- PLDDT Scores and Pocket IDs Table -->
      <div v-if="tableData.length > 0" class="q-mt-md">
        <div>Average PLDDT: {{ averagePLDDT.toFixed(2) }}</div>
        <h2>Confidence:</h2>
        <q-table
          :rows="tableData"
          :columns="columns"
          row-key="id"
        />
      </div>
    </q-page>
  </div>
</template>

<script lang="ts">
import { defineComponent } from 'vue';
import { QPage, QTable } from 'quasar';
import PdbViewer from 'src/components/PdbViewer.vue'; // Adjust the import path as necessary
import { useDrugDiscoveryStore } from 'src/features/drug_discovery/storage'; // Adjust the import path as necessary

export default defineComponent({
  name: 'UmolResult',
  components: {
    PdbViewer,
    QPage,
    QTable,
  },
  props: {
    experimentId: String,
    targetId: String,
    ligandId: String,
    jobId: String,
  },
  data() {
    return {
      predictedPdb: null as File | null,
      predictedSDF: null as File | null,
      plddtArray: [] as number[],
      tableData: [] as any[],
      columns: [
        { name: 'pocketIds', required: true, label: 'Pocket IDs', align: 'left', field: 'pocketIds', sortable: true },
        { name: 'plddt', label: 'PLDDT', align: 'left', field: 'plddt', sortable: true },
      ],
    };
  },
  computed: {
    averagePLDDT(): number {
      const total = this.plddtArray.reduce((acc, score) => acc + score, 0);
      return this.plddtArray.length ? total / this.plddtArray.length : 0;
    }
  },
  mounted() {
    this.loadUmolResults();
  },
  methods: {
    async loadUmolResults() {
      const store = useDrugDiscoveryStore();
      const resultData = await store.getUmolDockingJobResultData(this.experimentId!, this.targetId!, this.ligandId!, this.jobId!);
      if (resultData) {
        this.predictedPdb = new File([new Blob([resultData.predicted_pdb])], "predicted.pdb");
        this.predictedSDF = new File([new Blob([resultData.predicted_sdf])], "predicted.sdf");
        this.plddtArray = resultData.plddt_array;
      }

      const pocketData = await store.getJobPocketIds(this.experimentId!, this.targetId!, this.ligandId!, this.jobId!);
      if (pocketData && pocketData.pocket_ids) {
        // Assuming each pocket ID corresponds to a PLDDT score; adjust logic as needed
        this.tableData = pocketData.pocket_ids.map((id, index) => ({
          pocketIds: id,
          plddt: this.plddtArray[index] ?? null, // Handle cases where there might not be a corresponding PLDDT score
        }));
      }
    },
  },
});
</script>

<style scoped>
.q-table .q-table__container .q-table__middle {
  max-height: 300px; /* Adjust based on your preference */
  overflow-y: auto;
}
</style>
