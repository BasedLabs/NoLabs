<template>
  <q-card flat bordered v-if="!target.data">
    <q-card-section>
      <div class="text-caption">Protein sequence</div>
    </q-card-section>
    <q-separator />

    <q-card-section>
      <q-skeleton :type="text" />
    </q-card-section>
    <q-card-section>
      <div class="text-caption">3D Structure</div>
    </q-card-section>
    <q-skeleton height="200px" square />
    <q-separator />
  </q-card>
  <q-card v-if="target.data">
    <q-card-section>
      <div class="text-h6 q-pa-sm">Target: {{ this.target.metaData.target_name }}
        <q-btn round @click="changeTargetName"
                color="positive" size="sm" flat icon="edit"/></div>
      <q-card-section class="rounded-borders bg-black">
        <div class="text-h7 q-pl-md">Sequence: </div>
        <div class="fasta-sequence q-gutter-sm q-pa-md">
          <div v-for="row in wrappedFastaSequence" :key="row.index" class="row">
            <span
              v-for="(acid, index) in row.sequence"
              :key="index"
              :class="{ 'selected': isResidueSelected(row.startIndex + index) }"
              @click="toggleResidueSelection(row.startIndex + index)"
            >
              {{ acid }}
              <q-tooltip>{{`Residue ID: ${row.startIndex + index}` }}</q-tooltip>
            </span>
          </div>
        </div>
      </q-card-section>
      <div v-if="predictingFolding">
        <q-spinner color="info" size="lg" label="Loading 3D View..." />
      </div>
      <div class="q-pa-md">
        <q-select
          v-if="!this.target.data.pdbContents"
          v-model="selectedFoldingOption"
          :options="foldingOptions"
          label="Select folding method"
          outlined
          dense
          class="text-white"
        emit-value
        map-options
        ></q-select>
      </div>
      <q-btn v-if="!hasPdb" color="info" @click="predictStructure">
        Predict 3D structure
      </q-btn>
      <PdbViewer v-if="hasPdb && this.target.data.pdb_file" :pdb-file="this.target.data.pdb_file" :pocket-ids="this.target.data.pocketIds" :key="this.target.data.pocketIds" />
      <q-btn v-if="hasPdb && !bindingPocketAvailable" color="secondary" @click="fetchAndShowPocket">
        Show Binding Pocket
      </q-btn>
      <q-btn v-if="hasPdb && bindingPocketAvailable && !selectionMode" color="info" @click="toggleSelectionMode">
        Edit Binding Pocket
      </q-btn>
      <q-btn v-if="hasPdb && selectionMode" color="negative" @click="toggleSelectionMode">
        Cancel
      </q-btn>
      <q-btn
        v-if="selectionMode"
        color="positive"
        @click="sendSelectedResidues"
      >
        Confirm Selection
      </q-btn>
      <q-btn v-if="bindingPocketAvailable && !selectionMode" color="primary" @click="clearBindingPocket">
        Hide Binding Pocket
      </q-btn>

      <q-dialog v-model="showNoPocketModal">
        <q-card>
          <q-card-section>
            <div class="text-h6">Binding Pocket Not Found</div>
            <p>You don't have a binding pocket for this protein. Do you wish it to be predicted?</p>
            <div v-if="predictingPocket" class="q-mt-md">
              <q-spinner color="info" size="md" /> Predicting pocket...
            </div>
          </q-card-section>

          <q-card-actions align="right" v-if="!predictingPocket">
            <q-btn flat label="Close" color="white" @click="showNoPocketModal = false" />
            <q-btn flat label="Predict" color="white" @click="predictBindingPocket" />
          </q-card-actions>
        </q-card>
      </q-dialog>
    </q-card-section>
  </q-card>
</template>

<script>
import {QCard, QCardSection, QBtn, QDialog, QCardActions, Notify, QSpinnerOrbit, useQuasar} from "quasar";
import { useDrugDiscoveryStore } from "src/features/drug_discovery/storage";
import PdbViewer from "src/components/PdbViewer.vue";
import {obtainErrorResponse} from "../../../../api/errorWrapper";

export default {
  name: "TargetDetail",
  components: {
    QCard,
    QCardSection,
    QBtn,
    QDialog,
    QCardActions,
    PdbViewer
  },
  props: {
    experimentId: {
      type: String,
      required: true,
    },
    originalTarget: {
      type: Object,
      required: true,
    },
  },
  data() {
    return {
      selectedResidues: new Set(),
      showNoPocketModal: false,
      selectionMode: false,
      viewer: null,
      loading3DView: true,
      predictingFolding: false,
      predictingPocket: false,
      target: {
        metaData: {
          target_id: this.originalTarget.target_id,
          target_name: this.originalTarget.target_name
        },
        data: this.originalTarget.data
      },
      selectedFoldingOption: '',
      foldingOptions: [
        { label: 'Esmfold_light (up to 400 amino acids)', value: 'light' },
        { label: 'Esmfold', value: 'full' },
      ],
    };
  },
  methods: {
    async predictStructure() {
      const store = useDrugDiscoveryStore();
      if (this.target) {
        this.predictingFolding = true; // Start loading
        const predictionMethod = this.selectedFoldingOption === 'light' ? store.predictLightFoldingForTarget : store.predictFoldingForTarget;

        predictionMethod(this.experimentId, this.target.metaData.target_id)
          .then((response) => {
            if (response) {
              this.target.data.pdbContents = response;
              this.target.data.pdb_file = new File([new Blob([this.target.data.pdbContents])], "protein.pdb");
            }
          })
          .catch((error) => {
            console.error("Error predicting structure:", error);
          })
          .finally(() => {
            this.predictingFolding = false; // End loading
          });
      }
    },
    toggleSelectionMode() {
      this.selectionMode = !this.selectionMode;
    },
    toggleResidueSelection(index) {
      if (this.selectionMode) {
        let newPocketIds = [...this.target.data.pocketIds]; // Create a copy of the pocketIds array
        const pocketIndex = newPocketIds.indexOf(index);

        if (pocketIndex > -1) {
          newPocketIds.splice(pocketIndex, 1); // Remove the index
        } else {
          newPocketIds.push(index); // Add the index
        }

        this.target.data.pocketIds = newPocketIds; // Assign the new array back to pocketIds to trigger reactivity
      }
    },
    isResidueSelected(index) {
      return this.target.data.pocketIds && this.target.data.pocketIds.includes(index);
    },
    clearBindingPocket() {
      this.target.data.pocketIds = [];
    },
    async fetchAndShowPocket() {
      const store = useDrugDiscoveryStore();
      try {
        const pocketIds = await store.fetchPocketForTarget(this.experimentId, this.target.metaData.target_id);
        if (pocketIds && pocketIds.length > 0) {
          this.target.data.pocketIds = pocketIds;
        } else {
          this.showNoPocketModal = true; // Show modal to ask for prediction
        }
      } catch (error) {
        console.error('Error fetching binding pocket:', error);
      }
    },
    async predictBindingPocket() {
      const store = useDrugDiscoveryStore();
      this.predictingPocket = true; // Start the spinner
      try {
        const pocketIds = await store.predictPocketForTarget(this.experimentId, this.target.metaData.target_id);
        if (pocketIds) {
          this.target.data.pocketIds = pocketIds;
          this.bindingPocketAvailable = true;
        }
      } catch (error) {
        console.error('Error predicting binding pocket:', error);
      } finally {
        this.predictingPocket = false; // Stop the spinner
        this.showNoPocketModal = false; // Close the modal after prediction
      }
    },
    sendSelectedResidues() {
      const store = useDrugDiscoveryStore();
      if (this.target.data.pocketIds) {
        store
          .setPocketForTarget(
            this.experimentId,
            this.target.metaData.target_id,
            Array.from(this.target.data.pocketIds)
          )
          .then(() => {
            this.selectionMode = !this.selectionMode
          })
          .catch((error) => {
            console.error("Error setting binding pocket:", error);
          });
      }
    },
    changeTargetName() {
      const store = useDrugDiscoveryStore();
      this.quasar.dialog({
        color: 'positive',
        title: 'Prompt',
        message: 'Enter new target name',
        prompt: {
          model: this.target.metaData.target_name,
          required: true,
          type: 'text' // optional
        },
        cancel: true,
        persistent: true
      }).onOk(async data => {
        if (!data)
          return;
        this.quasar.loading.show({
          spinner: QSpinnerOrbit,
          message: 'Changing target name'
        });
        await store.changeTargetName(this.experimentId, this.target.metaData.target_id, data);
        this.target.metaData.target_name = data;
        this.quasar.loading.hide();
      });
    }
  },
  computed: {
    hasPdb() {
      return this.target.data.pdbContents != null;
    },
    bindingPocketAvailable() {
      return this.target.data.pocketIds && this.target.data.pocketIds.length > 0;
    },
    wrappedFastaSequence() {
      const sequence = this.target.data.proteinSequence || '';
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
  },
  mounted() {
    this.target.data.pdb_file = new File([new Blob([this.target.data.pdbContents])], "protein.pdb");
    this.quasar = useQuasar();
  },
  beforeUnmount() {
    if (this.viewer) {
      this.viewer.dispose();
    }
  },
};
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
