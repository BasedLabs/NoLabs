<template>
  <q-card flat bordered v-if="!protein.data">
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
  <q-card v-if="protein.data">
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
      <PdbViewer v-if="hasPdb && pdbFile" :pdb-file="pdbFile" :key="pdbFile.size" />
    </q-card-section>
  </q-card>
</template>

<script lang="ts">
import { QBtn, QCard, QCardSection, QSpinnerOrbit, QVueGlobals } from "quasar";
import PdbViewer from "src/components/PdbViewer.vue";
import { ProteinContentResponse } from "src/refinedApi/client";
import { defineComponent, PropType } from "vue";
import {useWorkflowStore} from "src/features/workflow/components/storage";

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
    originalProtein: {
      type: Object as PropType<ProteinContentResponse>,
      required: true,
    },
  },
  data() {
    return {
      loading3DView: true,
      protein: {
        metaData: {
          id: this.originalProtein.id,
          name: this.originalProtein.name,
        },
        data: {
          sequence: this.originalProtein.fasta_content,
          pdbContent: this.originalProtein.pdb_content
        }
      },
      quasar: null as QVueGlobals | null,
    };
  },
  methods: {
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
            console.error('Error changin protein name:', error);
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
