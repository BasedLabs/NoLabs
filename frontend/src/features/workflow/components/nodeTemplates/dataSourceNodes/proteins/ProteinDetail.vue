<template>
  <q-card flat bordered v-if="!protein.dataLoaded">
    <q-card-section>
      <div class="text-caption">Protein sequence</div>
    </q-card-section>
    <q-separator />

    <q-card-section>
      <q-skeleton type="text" />
    </q-card-section>
    <q-card-section>
      <div class="text-caption">3D Structure</div>
    </q-card-section>
    <q-skeleton height="200px" square />
    <q-separator />
  </q-card>
  <q-card v-if="protein.dataLoaded">
    <q-card-section>
      <div class="text-h6 q-pa-sm">Protein: {{ protein.metaData.name }}
        <q-btn round @click="changeProteinName" color="info" size="sm" flat icon="edit" />
      </div>
      <q-card-section class="rounded-borders bg-black">
        <div class="text-h7 q-pl-md">Sequence:</div>
        <div class="fasta-sequence q-gutter-sm q-pa-md">
          <div v-for="row in wrappedFastaSequence" :key="row.sequence" class="row">
            <span v-for="(acid, index) in row.sequence" :key="index">
              {{ acid }}
              <q-tooltip>{{ `Residue ID: ${row.startIndex + index}` }}</q-tooltip>
            </span>
          </div>
        </div>
      </q-card-section>
      <q-card-section v-if="protein.data.link">
        <q-btn flat color="primary" :href="protein.data.link" target="_blank">View More Info</q-btn>
      </q-card-section>
      <PdbViewer v-if="hasPdb && pdbFile" :pdb-file="pdbFile" :key="pdbFile.size" />
    </q-card-section>
  </q-card>
</template>

<script lang="ts">
import { QBtn, QCard, QCardSection, QSpinnerOrbit, QVueGlobals } from "quasar";
import PdbViewer from "src/components/PdbViewer.vue";
import { ProteinContentResponse } from "src/refinedApi/client";
import { defineComponent } from "vue";
import { getProteinContent } from "src/features/drug_discovery/refinedApi";
import { useWorkflowStore } from "src/features/drug_discovery/components/workflow/storage";

export default defineComponent({
  name: "ProteinDetail",
  components: {
    QCard,
    QCardSection,
    QBtn,
    PdbViewer
  },
  props: {
    experimentId: {
      type: String,
      required: true,
    },
    proteinId: {
      type: String,
      required: true,
    }
  },
  data() {
    return {
      loading3DView: true,
      protein: {
        metaData: {
          id: this.proteinId,
          name: '',
        },
        data: {
          sequence: '',
          pdbContent: '',
          link: ''
        },
        dataLoaded: false,
      },
      quasar: null as QVueGlobals | null,
    };
  },
  mounted() {
    this.loadProteinData();
  },
  methods: {
    async loadProteinData() {
      this.$q.loading.show({
        spinner: QSpinnerOrbit,
        message: 'Loading protein data'
      });

      try {
        const proteinContent: ProteinContentResponse = await getProteinContent(this.proteinId) as ProteinContentResponse;
        this.protein.metaData.name = proteinContent.name;
        this.protein.data.sequence = proteinContent.fasta_content || '';
        this.protein.data.pdbContent = proteinContent.pdb_content || '';
        this.protein.data.link = proteinContent.link || '';
        this.protein.dataLoaded = true;
      } catch (error) {
        console.error('Error loading protein data:', error);
      } finally {
        this.$q.loading.hide();
      }
    },
    changeProteinName() {
      this.$q.dialog({
        color: 'info',
        title: 'Prompt',
        message: 'Enter new protein name',
        prompt: {
          model: this.protein.metaData.name,
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
          message: 'Changing protein name'
        });
        try {
          await this.workflowStore.changeProteinName(this.protein.metaData.id, data);
          this.protein.metaData.name = data;
        } catch (error) {
          console.error('Error changing protein name:', error);
        }

        this.$q.loading.hide();
      });
    }
  },
  computed: {
    hasPdb() {
      return this.protein.data.pdbContent != null;
    },
    workflowStore() {
      return useWorkflowStore();
    },
    wrappedFastaSequence() {
      const sequence = this.protein.data.sequence || '';
      const wrapLength = 400; // Adjust based on your layout
      const wrapped = [];
      for (let i = 0; i < sequence.length; i += wrapLength) {
        wrapped.push({
          sequence: sequence.substring(i, i + wrapLength),
          startIndex: i,
        });
      }
      return wrapped;
    },
    pdbFile(): File {
      if (this.protein.data.pdbContent && this.protein.data.pdbContent.length > 0) {
        return new File([new Blob([this.protein.data.pdbContent])], this.protein.metaData.name + ".pdb");
      }

      return new File([], "empty_protein.pdb");
    }
  }
})
</script>

<style>
.fasta-sequence {
  display: flex;
  flex-wrap: wrap;
}

.fasta-sequence .row {
  display: flex;
  flex-wrap: wrap;
}

.fasta-sequence span {
  padding: 2px;
  margin-right: 1px;
  cursor: pointer;
}

.fasta-sequence span.selected {
  background-color: blue;
  color: white;
}
</style>
