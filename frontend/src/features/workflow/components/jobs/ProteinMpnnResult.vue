<template>
  <div class="protein-mpnn-result">
    <div class="generated-sequences-table">
      <h3>Generated Sequences</h3>
      <q-table :rows="sequences" :columns="sequenceColumns" row-key="id">
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

interface SequenceEntry {
  id: number;
  sequence: string;
  score?: number;
  global_score?: number;
  T?: number;
  sample?: number;
  seq_recovery?: number;
  [key: string]: number | string | undefined; // For any additional parameters
}

export default defineComponent({
  name: 'ProteinMPNNResult',
  props: {
    job: Object as () => { fasta_contents: string[] },
    protein: Object,
  },
  data() {
    return {
      sequences: [] as SequenceEntry[],
      sequenceColumns: [
        { name: 'id', label: 'ID', field: 'id', align: 'left' },
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
    this.parseFastaContents();
  },
  methods: {
    parseFastaContents() {
      const fastaContents = this.job.fasta_contents;
      const sequences = [];

      fastaContents.forEach((content) => {
        // Split the content by '>'
        const entries = content.split('>').filter(entry => entry.trim() !== '');
        entries.forEach((entry, index) => {
          const lines = entry.split('\n').filter(line => line.trim() !== '');
          const header = lines[0];
          const sequence = lines.slice(1).join('').replace(/\s+/g, '');

          // Parse parameters from header
          const params = this.parseHeader(header);

          sequences.push({
            id: index + 1,
            sequence,
            ...params,
          });
        });
      });

      this.sequences = sequences;
    },
    parseHeader(header: string) {
      const params: { [key: string]: number | string } = {};

      // Use regular expressions to extract key-value pairs
      const regex = /(\w+)=([^\s,]+)/g;
      let match;
      while ((match = regex.exec(header)) !== null) {
        params[match[1]] = match[2].replace(/^\[|\]$/g, ''); // Remove brackets if present
      }

      // Additionally, extract 'T', 'sample', etc., which might be at the start without keys
      const additionalParamsRegex = /^(T=[^\s,]+|sample=[^\s,]+)/g;
      const additionalParams = header.match(additionalParamsRegex);
      if (additionalParams) {
        additionalParams.forEach(param => {
          const [key, value] = param.split('=');
          params[key] = value;
        });
      }

      return params;
    },
    downloadFasta() {
      const fastaContent = this.job.fasta_contents.join('\n');
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
