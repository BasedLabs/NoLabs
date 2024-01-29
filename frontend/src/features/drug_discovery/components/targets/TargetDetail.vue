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
      <div class="text-h6">Target: {{ this.target.metaData.target_name }}</div>
      <div class="amino-acid-sequence">
        <span
          v-for="(acid, index) in fastaSequence"
          :key="index"
          :class="{ selected: selectedResidues.has(index) }"
          @click="toggleFastaResidueSelection(index)"
        >
          {{ acid }}
        </span>
      </div>
      <div v-if="predictingFolding" class="loading-indicator">
        <q-spinner color="primary" label="Loading 3D View..." />
      </div>
      <q-btn v-if="!hasPdb" color="accent" @click="predictStructure">
        Predict 3D structure
      </q-btn>
      <div v-else ref="viewerContainer" class="protein-viewer" />
      <q-btn v-if="hasPdb && !bindingPocketAvailable" color="secondary" @click="fetchAndShowPocket">
        Show Binding Pocket
      </q-btn>
      <q-btn v-if="hasPdb" color="accent" @click="toggleSelectionMode">
        Edit Binding Pocket
      </q-btn>
      <q-btn
        v-if="selectionMode"
        color="positive"
        @click="sendSelectedResidues"
      >
        Confirm Selection
      </q-btn>
      <q-dialog v-model="showNoPocketModal">
        <q-card>
          <q-card-section>
            <div class="text-h6">Binding Pocket Not Found</div>
            <p>You don't have a binding pocket for this protein. Do you wish it to be predicted?</p>
          </q-card-section>

          <q-card-actions align="right">
            <q-btn flat label="Close" color="primary" @click="showNoPocketModal = false" />
            <q-btn flat label="Predict" color="primary" @click="predictBindingPocket" />
          </q-card-actions>
        </q-card>
      </q-dialog>
    </q-card-section>
  </q-card>
</template>

<script>
import { QCard, QCardSection, QBtn, QDialog, QCardActions } from "quasar";
import { useDrugDiscoveryStore } from "src/features/drug_discovery/storage";
import * as NGL from "ngl";

export default {
  name: "TargetDetail",
  components: {
    QCard,
    QCardSection,
    QBtn,
    QDialog,
    QCardActions,
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
      target: {
        metaData: {
          target_id: this.originalTarget.target_id,
          target_name: this.originalTarget.target_name
        },
        data: this.originalTarget.data
      },
    };
  },
  methods: {
    createViewer() {
      if (this.viewer || !this.target.data.pdbContents) {
        this.loading3DView = false; // Ensure loading state is cleared even if conditions aren't met
        return;
      }
      setTimeout(async () => {
        this.viewer = new NGL.Stage(this.$refs.viewerContainer);

        const stringBlob = new Blob([this.target.data.pdbContents], { type: "text/plain" });
        const proteinFile = new File([stringBlob], 'protein.pdb', { type: 'text/plain' });
        this.viewer.loadFile(proteinFile, { defaultRepresentation: true }).then((component) => {
          component.autoView();
          this.loading3DView = false;
        });
      }, 100);
    },
    predictStructure() {
      const store = useDrugDiscoveryStore();
      if (this.target) {
        this.predictingFolding = true; // Start loading
        store
          .predictFoldingForTarget(this.experimentId, this.target.metaData.target_id)
          .then((response) => {
            this.target.data.pdbContents = response;
            // Re-render 3D structure with new PDB data, if applicable
            this.createViewer();
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
    toggleFastaResidueSelection(index) {
      if (this.selectedResidues.has(index)) {
        this.selectedResidues.delete(index);
      } else {
        this.selectedResidues.add(index);
      }
      this.highlightSelectedResidues(); // Update the 3D view to reflect the selection
    },
    async fetchAndShowPocket() {
      const store = useDrugDiscoveryStore();
      try {
        const pocketIds = await store.fetchPocketForTarget(this.experimentId, this.target.metaData.target_id);
        if (pocketIds && pocketIds.length > 0) {
          this.target.data.pocketIds = pocketIds;
          this.highlightSelectedResidues();
        } else {
          this.showNoPocketModal = true; // Show modal to ask for prediction
        }
      } catch (error) {
        console.error('Error fetching binding pocket:', error);
      }
    },
    async predictBindingPocket() {
      const store = useDrugDiscoveryStore();
      this.target.data.pocketIds = await store.predictPocketForTarget(this.experimentId,
        this.target.metaData.target_id);
      this.bindingPocketAvailable = true;
      this.highlightSelectedResidues(); // Update 3D view to reflect the new selection
    },
    highlightSelectedResidues() {
      if (!this.viewer || !this.target.data.pdbContents) {
        return;
      }
      // Update the 3D viewer to highlight the selected residues
      const component = this.viewer.compList[0];
      component.removeAllRepresentations();
      component.addRepresentation("cartoon");

      // Convert pocketIds to a format suitable for NGL selection (1-indexed)
      const selectionString = this.target.data.pocketIds
        .map((id) => (id + 1).toString())
        .join(" or ");
      component.addRepresentation("ball+stick", {
        sele: selectionString,
        color: "blue",
      });
      component.autoView();
    },
    sendSelectedResidues() {
      const store = useDrugDiscoveryStore();
      if (this.target.data.pocketIds) {
        store
          .setBindingPocket(
            this.target.metaData.target_id,
            Array.from(this.selectedResidues)
          )
          .then(() => {
            // Handle success
          })
          .catch((error) => {
            console.error("Error setting binding pocket:", error);
          });
      }
    },
  },
  computed: {
    hasPdb() {
      return this.target.data.pdbContents != null;
    },
    bindingPocketAvailable() {
      return this.target.data.pocketIds && this.target.data.pocketIds.length > 0;
    },
    fastaSequence() {
      return this.target.data.proteinSequence ? this.target.data.proteinSequence.split('') : [];
    },
  },
  mounted() {
   this.createViewer();
  },
  beforeUnmount() {
    if (this.viewer) {
      this.viewer.dispose();
    }
  },
};
</script>

<style>
.amino-acid-sequence {
  overflow-x: auto;
  /* Allows horizontal scrolling for long sequences */
  white-space: nowrap;
  /* Ensures the sequence is in one line */
  font-family: monospace;
  /* Gives a uniform look to the sequence characters */
  margin-bottom: 10px;
  /* Spacing below the sequence */
  border: 1px solid #ddd;
  /* Border around the sequence */
  border-radius: 5px;
  /* Rounded corners for the border */
  padding: 10px;
  /* Padding inside the border */
  background-color: #000000;
  /* Light background for the sequence */
}

.amino-acid-sequence span {
  display: inline-block;
  /* Each amino acid is an inline-block element */
  padding: 2px 5px;
  /* Padding around each amino acid */
  cursor: pointer;
  /* Indicates the amino acid can be clicked */
  border-radius: 3px;
  /* Rounded corners for each amino acid span */
  margin-right: 2px;
  /* Space between amino acids */
}

.amino-acid-sequence span.selected {
  background-color: blue;
  /* Background color for selected amino acids */
  color: white;
  /* Text color for selected amino acids */
}

.loading-indicator {
  display: flex;
  align-items: center;
  justify-content: center;
  padding: 10px;
}


.protein-viewer {
  width: 100%;
  height: 400px;
}
</style>
