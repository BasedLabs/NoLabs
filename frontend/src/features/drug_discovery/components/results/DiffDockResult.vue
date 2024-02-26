<template>
  <div class="diff-dock-result">
    <div class="predicted-ligands-table">
      <h3>Predicted Ligands</h3>
      <q-table :rows="predictedLigands" :columns="ligandColumns" row-key="predicted_ligand_file_name">
        <template v-slot:body="props">
          <q-tr
            :props="props"
            @click="() => handleRowClick(props.row)"
            :class="{ 'selected-row': props.row.predicted_ligand_file_name === selectedLigandFileName }"
          >
            <q-td v-for="col in ligandColumns" :key="col.name" :props="props">
              {{ props.row[col.field] }}
            </q-td>
          </q-tr>
        </template>
      </q-table>
    </div>
    <div class="pdb-viewer">
      <h3>Predicted Structure</h3>
      <PdbViewer v-if="predictedPdb" :pdb-file="predictedPdb" :sdf-file="predictedSDF" :key="predictedSDF" />
    </div>
  </div>
</template>

<script lang="ts">
import { defineComponent } from 'vue';
import PdbViewer from 'src/components/PdbViewer.vue';
import {DiffDockLigandMetaData} from "src/api/client/models/DiffDockLigandMetaData.ts";
import {useDrugDiscoveryStore} from 'src/features/drug_discovery/storage';

export default defineComponent({
  name: 'DiffDockResult',
  components: { PdbViewer },
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
      selectedLigandFileName: null as string | null,
      predictedLigands: [] as DiffDockLigandMetaData[],
      ligandColumns: [
        { name: 'predicted_ligand_file_name', label: 'Ligand File Name', field: 'predicted_ligand_file_name', align: 'left' },
        {
          name: 'confidence',
          label: 'Confidence',
          field: 'confidence',
          align: 'left',
          format: (val: number | null) => val?.toFixed(2) ?? 'N/A' // Specify type for val
        },
        {
          name: 'scored_affinity',
          label: 'Scored Affinity',
          field: 'scored_affinity',
          align: 'left',
          format: (val: number) => val.toFixed(2) // Specify type for val
        },
        {
          name: 'minimized_affinity',
          label: 'Minimized Affinity',
          field: 'minimized_affinity',
          align: 'left',
          format: (val: number) => val.toFixed(2) // Specify type for val
        },
      ],
    };
  },
  async mounted() {
    await this.loadDiffDockResults();

    if (this.predictedLigands.length > 0) {
      this.loadLigandSDF(this.predictedLigands[0].predicted_ligand_file_name);
    }
  },
  methods: {
    async loadDiffDockResults() {
      const store = useDrugDiscoveryStore();
      const result = await store.getDiffDockDockingJobResultData(this.experimentId!, this.targetId!, this.ligandId!, this.jobId!);
      if (result) {
        this.predictedPdb = new File([new Blob([result.predicted_pdb])], "protein.pdb");
        this.predictedLigands = result.predicted_ligands;
      }
    },
    async loadLigandSDF(predictedLigandFileName: string) {
      const store = useDrugDiscoveryStore();
      const sdfResp = await store.getDiffDockLigandSdf(
        this.experimentId!,
        this.targetId!,
        this.ligandId!,
        this.jobId!,
        predictedLigandFileName
      );
      if (sdfResp?.sdf_contents) {
        this.predictedSDF = new File([new Blob([sdfResp.sdf_contents])], predictedLigandFileName);
        this.selectedLigandFileName = predictedLigandFileName;
      }
    },
    handleRowClick(row: DiffDockLigandMetaData) {
      this.loadLigandSDF(row.predicted_ligand_file_name);
    },
  },
});
</script>

<style scoped>
.selected-row {
  background-color: #1976d2; /* Adjust the color as needed */
}
</style>
