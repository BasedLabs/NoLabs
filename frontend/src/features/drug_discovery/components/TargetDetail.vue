<template>
    <q-page class="target-detail-page">
        <q-card v-if="target">
            <q-card-section>
                <div class="amino-acid-sequence">
                    <span v-for="(acid, index) in fastaSequence" :key="index"
                        :class="{ 'selected': selectedResidues.has(index) }" @click="toggleFastaResidueSelection(index)">
                        {{ acid }}
                    </span>
                </div>

                <div v-if="hasPdb">
                    <!-- 3D Protein Structure View Here -->
                    <!-- Placeholder for 3D view component or element -->
                </div>
                <q-btn v-else color="primary" @click="predictStructure">Predict 3D structure</q-btn>

                <q-btn v-if="hasPdb" color="secondary" @click="toggleDisplayBindingPocket">
                    {{ bindingPocketAvailable ? 'Display Binding Pocket' : 'Predict Binding Pocket' }}
                </q-btn>

                <!-- Modal for Binding Pocket Prediction -->
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
                <q-btn v-if="hasPdb" color="accent" @click="toggleSelectionMode">
                    Edit Binding Pocket
                </q-btn>

                <q-btn v-if="selectionMode" color="positive" @click="sendSelectedResidues">
                    Confirm Selection
                </q-btn>
            </q-card-section>
        </q-card>
    </q-page>
</template>
<script>
import { defineComponent, computed, ref } from 'vue';
import { QCard, QCardSection, QBtn, QDialog, QCardActions } from 'quasar';

export default defineComponent({
    name: 'TargetDetail',
    components: {
        QCard,
        QCardSection,
        QBtn,
        QDialog,
        QCardActions,
    },
    props: {
        target: Object, // Assuming target is an object containing necessary data
    },
    setup(props) {
        const selectedResidues = ref(new Set());
        const showNoPocketModal = ref(false);
        const selectionMode = ref(false);

        const fastaSequence = computed(() => props.target?.fasta.split('') || []);

        const hasPdb = computed(() => !!props.target?.pdb);
        const bindingPocketAvailable = computed(() => !!props.target?.bindingPocket);

        const toggleFastaResidueSelection = (index) => {
            if (selectedResidues.value.has(index)) {
                selectedResidues.value.delete(index);
            } else {
                selectedResidues.value.add(index);
            }
        };

        const predictStructure = () => {
            // Logic for predicting the structure
        };

        const toggleDisplayBindingPocket = () => {
            if (bindingPocketAvailable.value) {
                // Logic to display the binding pocket
            } else {
                show
                NoPocketModal.value = true; // Show modal to predict the binding pocket
            }
        };
        const predictBindingPocket = () => {
            // Logic to predict the binding pocket
            showNoPocketModal.value = false; // Close the modal after initiating prediction
        };

        const toggleSelectionMode = () => {
            selectionMode.value = !selectionMode.value;
            // Add additional logic for entering/exiting selection mode if needed
        };

        const sendSelectedResidues = () => {
            // Logic to handle sending selected residues for binding pocket
        };

        // Additional methods and computed properties as needed

        return {
            fastaSequence,
            hasPdb,
            bindingPocketAvailable,
            selectedResidues,
            showNoPocketModal,
            selectionMode,
            toggleFastaResidueSelection,
            predictStructure,
            toggleDisplayBindingPocket,
            predictBindingPocket,
            toggleSelectionMode,
            sendSelectedResidues,
        };
    },
});
</script>
  
<style>
.target-detail-page {
    /* Adjust layout and spacing for the target detail page */
}

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
    background-color: #f8f8f8;
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
}</style>