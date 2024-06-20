<template>
    <div class="diff-dock-result">
      <div class="predicted-ligands-table">
        <h3>Predicted Ligands</h3>
        <q-table :rows="predictedLigands" :columns="ligandColumns" row-key="complex_id">
          <template v-slot:body="props">
            <q-tr
              :props="props"
              @click="handleRowClick(props.row)"
              :class="{ 'selected-row': props.row.complex_id === selectedLigandFileName }"
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
  import {
    nolabs__application__use_cases__diffdock__api_models__JobResult,
    ProteinContentResponse
  } from "../../../../refinedApi/client";

  export default defineComponent({
    name: 'DiffDockResult',
    components: { PdbViewer },
    props: {
      job: Object as () => { result: Array<nolabs__application__use_cases__diffdock__api_models__JobResult> },
      protein: Object as () => ProteinContentResponse,
    },
    data() {
      return {
        predictedPdb: null as File | null,
        predictedSDF: null as File | null,
        selectedLigandFileName: null as string | null,
        predictedLigands: [] as nolabs__application__use_cases__diffdock__api_models__JobResult[],
        ligandColumns: [
          { name: 'complex_id', label: 'Complex ID', field: 'complex_id', align: 'left' },
          {
            name: 'confidence',
            label: 'Confidence',
            field: 'confidence',
            align: 'left',
            format: (val: number | null) => val?.toFixed(2) ?? 'N/A',
          },
          {
            name: 'scored_affinity',
            label: 'Scored Affinity',
            field: 'scored_affinity',
            align: 'left',
            format: (val: number) => val.toFixed(2),
          },
          {
            name: 'minimized_affinity',
            label: 'Minimized Affinity',
            field: 'minimized_affinity',
            align: 'left',
            format: (val: number) => val.toFixed(2),
          },
        ],
      };
    },
    async mounted() {
      this.predictedLigands = this.job.result;
      if (this.predictedLigands.length > 0) {
        this.loadLigandSDF(this.predictedLigands[0].complex_id);
      }
    },
    methods: {
      async loadLigandSDF(complexId: string) {
        const ligand = this.predictedLigands.find(ligand => ligand.complex_id === complexId);
        if (ligand && this.protein) {
          this.predictedSDF = new File([new Blob([ligand.sdf_content])], complexId + ".sdf");
          this.predictedPdb = new File([new Blob([this.protein.pdb_content!!])], this.protein.name + ".pdb");
          this.selectedLigandFileName = complexId;
        }
      },
      handleRowClick(row: nolabs__application__use_cases__diffdock__api_models__JobResult) {
        this.loadLigandSDF(row.complex_id);
      },
    },
  });
  </script>

  <style scoped>
  .selected-row {
    background-color: #1976d2; /* Adjust the color as needed */
  }
  </style>
