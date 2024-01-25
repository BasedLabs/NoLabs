<template>
    <q-card flat bordered v-if="!target.data">
        <q-card-section>
            <div class="text-caption">Protein sequence</div>
        </q-card-section>
        <q-separator />

        <q-card-section>
            <q-skeleton :type="type" />
        </q-card-section>
        <q-card-section>
            <div class="text-caption">3D Structure</div>
        </q-card-section>
        <q-skeleton height="200px" square />
        <q-separator />
    </q-card>
    <q-card v-if="target.data">
        <q-card-section>
            <div class="text-h6">Target: {{ this.target.metaData.name }}</div>
            <div class="amino-acid-sequence">
                <span v-for="(acid, index) in fastaSequence" :key="index"
                    :class="{ 'selected': selectedResidues.has(index) }" @click="toggleFastaResidueSelection(index)">
                    {{ acid }}
                </span>
            </div>
            <div v-if="predictingFolding">
                <q-spinner color="primary" label="Predicting folding..." />
            </div>
            <div v-else>
                <div v-if="hasPdb" ref="viewerContainer" class="protein-viewer"></div>
                <q-btn v-else color="primary" @click="predictStructure">Predict 3D structure</q-btn>
                <q-btn v-if="hasPdb" color="secondary" @click="toggleDisplayBindingPocket">
                    {{ bindingPocketAvailable ? 'Display Binding Pocket' : 'Predict Binding Pocket' }}
                </q-btn>
                <q-btn v-if="hasPdb" color="accent" @click="toggleSelectionMode">
                    Edit Binding Pocket
                </q-btn>

                <q-btn v-if="selectionMode" color="positive" @click="sendSelectedResidues">
                    Confirm Selection
                </q-btn>
            </div>
            <q-dialog v-model="showNoPocketModal">
                <q-card>
                    <q-card-section>
                        <div class="text-h6">Binding Pocket Not Found</div>
                        <p>You don't have a binding pocket for this protein. Do you wish it to be predicted?</p>
                    </q-card-section>

                    <q-card-actions align="right">
                        <q-btn flat label="Close" color="primary" v-close-popup />
                        <q-btn flat label="Predict" color="primary" @click="predictBindingPocket" />
                    </q-card-actions>
                </q-card>
            </q-dialog>
        </q-card-section>
    </q-card>
</template>


<script>
import { QCard, QCardSection, QBtn, QDialog, QCardActions } from 'quasar';
import { useDrugDiscoveryStore } from '../../storage';
import * as NGL from 'ngl';

export default {
    name: 'TargetDetail',
    components: {
        QCard,
        QCardSection,
        QBtn,
        QDialog,
        QCardActions,
    },
    props: {
        targetMetadata: {
            type: Object,
            required: false,
        }
    },
    data() {
        return {
            selectedResidues: new Set(),
            showNoPocketModal: false,
            selectionMode: false,
            viewer: null,
            target: {
                metaData: this.targetMetadata,
                data: null
            }
        };
    },
    methods: {
        createViewer() {
            if (this.viewer || !this.target.pdbContents) {
                return;
            }

            // Initialize NGL viewer
            this.viewer = new NGL.Stage(this.$refs.viewerContainer);

            // Load PDB data into viewer
            const stringBlob = new Blob([this.target.pdbContents], { type: 'text/plain' });
            this.viewer.loadFile(stringBlob, { ext: 'pdb' }).then((component) => {
                component.addRepresentation('cartoon', { color: 'sstruc' });
                component.autoView();
            });
        },
        predictStructure() {
            const store = useDrugDiscoveryStore();
            if (this.target) {
                this.predictingFolding = true; // Start loading
                store.predictFoldingForTarget(this.experimentId, this.target.targetId).then(response => {
                    this.target.pdbContents = response.pdbContents;
                    // Re-render 3D structure with new PDB data, if applicable
                }).catch(error => {
                    console.error('Error predicting structure:', error);
                }).finally(() => {
                    this.predictingFolding = false; // End loading
                });
            }
        },
        toggleDisplayBindingPocket() {
            if (this.bindingPocketAvailable) {
                // Logic to display the binding pocket
            } else {
                this.showNoPocketModal = true;
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
        predictBindingPocket() {
            // Mock data for binding pocket prediction
            this.localTarget.pocketIds = Array.from({ length: 41 }, (_, i) => 50 + i); // [50, 51, ..., 90]
            this.highlightSelectedResidues(); // Update 3D view to reflect the new selection
        },
        highlightSelectedResidues() {
            if (!this.viewer || !this.localTarget.pdbContents) {
                return;
            }

            // Update the 3D viewer to highlight the selected residues
            const component = this.viewer.compList[0];
            component.removeAllRepresentations();
            component.addRepresentation('cartoon', { color: 'sstruc' });

            // Convert pocketIds to a format suitable for NGL selection (1-indexed)
            const selectionString = this.localTarget.pocketIds.map(id => (id + 1).toString()).join(' or ');
            component.addRepresentation('ball+stick', {
                sele: selectionString,
                color: 'blue',
            });

            component.autoView();
        },
        sendSelectedResidues() {
            const store = useDrugDiscoveryStore();
            if (this.localTarget) {
                store.setBindingPocket(this.localTarget.targetId, Array.from(this.selectedResidues)).then(() => {
                    // Handle success
                }).catch(error => {
                    console.error('Error setting binding pocket:', error);
                });
            }
        },
    },
    mounted() {
        // Any additional logic needed on mount, if applicable
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

.q-dialog .q-card {
    width: 400px;
    /* Width of the modal dialog */
}

.q-dialog .q-card-section {
    color: black;
    /* Text color inside the modal */
}

.q-dialog .q-card-actions {
    justify-content: flex-end;
    /* Align modal actions (buttons) to the right */
}

.protein-viewer {
    width: 100%;
    height: 400px;
}
</style>