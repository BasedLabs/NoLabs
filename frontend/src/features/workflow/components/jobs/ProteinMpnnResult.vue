<template>
  <div class="protein-mpnn-result">
    <div class="generated-sequences-table">
      <h3>Generated Sequences</h3>
      <q-table :rows="sequences" :columns="sequenceColumns" row-key="sequence_id">
        <!-- Custom body slot for better control -->
        <template v-slot:body="props">
          <q-tr :props="props">
            <q-td v-for="col in sequenceColumns" :key="col.name" :props="props">
              <div v-if="col.name === 'sequence'" class="sequence-cell">{{ props.row[col.field] }}</div>
              <div v-else>{{ props.row[col.field] }}</div>
            </q-td>
          </q-tr>
        </template>
      </q-table>
      <!-- Download Button -->
      <q-btn color="primary" @click="downloadFasta" label="Download FASTA"></q-btn>
    </div>
  </div>
</template>

<script lang="ts">
import { defineComponent } from 'vue';
import { QTable, QTr, QTd, QBtn } from 'quasar';
import {
  nolabs__application__proteinmpnn__api_models__JobResponse,
  nolabs__application__proteinmpnn__api_models__JobResult,
} from 'src/refinedApi/client';

export default defineComponent({
  name: 'ProteinMPNNResult',
  props: {
    job: {
      type: Object as () => nolabs__application__proteinmpnn__api_models__JobResponse,
      required: true,
    },
    protein: Object,
  },
  data() {
    return {
      sequences: [] as nolabs__application__proteinmpnn__api_models__JobResult[],
      sequenceColumns: [
        { name: 'sequence_id', label: 'ID', field: 'sequence_id', align: 'left' },
        { name: 'score', label: 'Score', field: 'score', align: 'left' },
        { name: 'global_score', label: 'Global Score', field: 'global_score', align: 'left' },
        { name: 'T', label: 'Temperature', field: 'T', align: 'left' },
        { name: 'sample', label: 'Sample', field: 'sample', align: 'left' },
        { name: 'seq_recovery', label: 'Seq Recovery', field: 'seq_recovery', align: 'left' },
        { name: 'sequence', label: 'Sequence', field: 'sequence', align: 'left' },
      ],
    };
  },
  mounted() {
    this.loadSequences();
  },
  methods: {
    loadSequences() {
      if (this.job && this.job.result) {
        this.sequences = this.job.result;
      }
    },
    downloadFasta() {
      if (!this.job || !this.job.result) return;

      const fastaContent = this.job.result
        .map((result) => result.fasta_content)
        .join('\n');

      const blob = new Blob([fastaContent], { type: 'text/plain;charset=utf-8' });
      const url = URL.createObjectURL(blob);

      const link = document.createElement('a');
      link.href = url;
      link.download = 'protein_mpnn_results.fasta';
      link.click();

      URL.revokeObjectURL(url);
    },
  },
  components: {
    QTable,
    QTr,
    QTd,
    QBtn,
  },
});
</script>

<style scoped>
.sequence-cell {
  max-width: 400px;
  word-wrap: break-word;
  white-space: pre-wrap;
}
</style>
