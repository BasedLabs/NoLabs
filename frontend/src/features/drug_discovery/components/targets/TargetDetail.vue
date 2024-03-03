<template>
  <q-card flat bordered v-if="!target.data">
    <q-card-section>
      <div class="text-caption">Protein sequence</div>
    </q-card-section>
    <q-separator/>

    <q-card-section>
      <q-skeleton type="text"/>
    </q-card-section>
    <q-card-section>
      <div class="text-caption">3D Structure</div>
    </q-card-section>
    <q-skeleton height="200px" square/>
    <q-separator/>
  </q-card>
  <q-card v-if="target.data">
    <q-card-section>
      <div class="text-h6 q-pa-sm">Target: {{ this.target.metaData.target_name }}
        <q-btn round @click="changeTargetName"
               color="info" size="sm" flat icon="edit"/>
      </div>
      <q-card-section class="rounded-borders bg-black">
        <div class="text-h7 q-pl-md">Sequence:</div>
        <div class="fasta-sequence q-gutter-sm q-pa-md">
          <div v-for="row in wrappedFastaSequence" :key="row.index" class="row">
            <span
                v-for="(acid, index) in row.sequence"
                :key="index"
                :class="{ 'selected': isResidueSelected(row.startIndex + index) }"
                @click="toggleResidueSelection(row.startIndex + index)"
            >
              {{ acid }}
              <q-tooltip>{{ `Residue ID: ${row.startIndex + index}` }}</q-tooltip>
            </span>
          </div>
        </div>
      </q-card-section>
      <div v-if="predictingFolding">
        <q-spinner color="info" size="lg" label="Loading 3D View..."/>
      </div>
      <div class="q-pa-md">
        <q-select
            v-if="!this.target.data.pdbContents"
            v-model="selectedFoldingOption"
            :options="foldingOptionsSelect"
            label="Select folding method"
            outlined
            dense
            class="text-white"
            emit-value
            map-options
        ></q-select>
      </div>
      <q-btn v-if="!hasPdb" color="primary" @click="predictStructure">
        Predict 3D structure
      </q-btn>
      <PdbViewer v-if="hasPdb && this.pdbFile" :pdb-file="this.pdbFile"
                 :pocket-ids="this.target.data.pocketIds" :key="this.target.data.pocketIds"/>
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
          color="info"
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
              <q-spinner color="info" size="md"/>
              Predicting pocket...
            </div>
          </q-card-section>

          <q-card-actions align="right" v-if="!predictingPocket">
            <q-btn flat label="Close" color="white" @click="showNoPocketModal = false"/>
            <q-btn flat label="Predict" color="white" @click="predictBindingPocket"/>
          </q-card-actions>
        </q-card>
      </q-dialog>
    </q-card-section>
  </q-card>
</template>

<script lang="ts">
import {Notify, QBtn, QCard, QCardActions, QCardSection, QDialog, QSpinnerOrbit, QVueGlobals} from "quasar";
import {useDrugDiscoveryStore} from "src/features/drug_discovery/storage";
import PdbViewer from "src/components/PdbViewer.vue";
import {defineComponent, PropType} from "vue";

export default defineComponent({
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
      type: Object as PropType<{
        target_id: string,
        target_name: string,
        data: {
          pocketIds: number[],
          pdbContents: string | undefined | null,
          proteinSequence: string | undefined | null
        }
      }>,
      required: true,
    },
  },
  data() {
    return {
      selectedResidues: new Set(),
      showNoPocketModal: false,
      selectionMode: false,
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
      foldingOptions: {
        esmfoldLight: 'light',
        esmfold: 'esmfold',
        rosettafold: 'rosettafold'
      },
      quasar: null as QVueGlobals | null,
    };
  },
  methods: {
    async predictStructure() {
      const store = useDrugDiscoveryStore();
      if (this.target) {
        let predictionMethod: null | ((experimentId: string, targetId: string) => Promise<string | null | undefined>) = null;

        if (this.selectedFoldingOption === this.foldingOptions.esmfoldLight) {
          predictionMethod = store.predictEsmLightFoldingForTarget;
        }

        if (this.selectedFoldingOption === this.foldingOptions.esmfold) {
          predictionMethod = store.predictEsmFoldingForTarget;
        }

        if (this.selectedFoldingOption === this.foldingOptions.rosettafold) {
          predictionMethod = store.predictRoseTTAFold;
        }

        if (predictionMethod === null) {
          Notify.create({
            type: "negative",
            closeBtn: 'Close',
            message: 'Select folding method'
          });
          return;
        }

        try {
          this.predictingFolding = true; // Start loading
          this.target.data.pdbContents = await predictionMethod(this.experimentId, this.target.metaData.target_id);
        } catch (error) {
          console.error("Error prediction structure", error);
        }

        this.predictingFolding = false;
      }
    },
    toggleSelectionMode() {
      this.selectionMode = !this.selectionMode;
    },
    toggleResidueSelection(index: number) {
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
    isResidueSelected(index: number) {
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
      }

      this.predictingPocket = false; // Stop the spinner
      this.showNoPocketModal = false; // Close the modal after prediction
    },
    async sendSelectedResidues() {
      const store = useDrugDiscoveryStore();
      if (this.target.data.pocketIds) {
        try {
          await store
              .setPocketForTarget(
                  this.experimentId,
                  this.target.metaData.target_id,
                  Array.from(this.target.data.pocketIds)
              );
          this.selectionMode = !this.selectionMode;
        } catch (error) {
          console.error("Error setting binding pocket:", error);
        }
      }
    },
    changeTargetName() {
      const store = useDrugDiscoveryStore();
      this.$q.dialog({
        color: 'info',
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
        this.$q.loading.show({
          spinner: QSpinnerOrbit,
          message: 'Changing target name'
        });
        await store.changeTargetName(this.experimentId, this.target.metaData.target_id, data);
        this.target.metaData.target_name = data;
        this.$q.loading.hide();
      });
    }
  },
  computed: {
    hasPdb() {
      return this.target.data.pdbContents != null;
    },
    foldingOptionsSelect() {
      return [
        {label: 'Esmfold_light (up to 400 amino acids)', value: this.foldingOptions.esmfoldLight},
        {label: 'Esmfold', value: this.foldingOptions.esmfold},
        {label: 'RoseTTAFold', value: this.foldingOptions.rosettafold},
      ];
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
    pdbFile(): File {
      if(this.target.data.pdbContents && this.target.data.pdbContents.length > 0){
        return new File([new Blob([this.target.data.pdbContents])], "protein.pdb");
      }

      return new File([], "protein.pdb");
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
