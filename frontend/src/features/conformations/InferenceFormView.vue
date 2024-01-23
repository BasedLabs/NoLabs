<script lang="ts">
import {defineComponent, PropType} from 'vue'

interface OnSubmitType {
  (contig: string, hotspots: string, numberOfDesigns: number, timesteps: number, pdbFile: File): Promise<void>;
}

interface DataType {
  contig: string;
  hotspots: string;
  numberOfDesigns: number;
  timesteps: number;
  pdbFile: File | null;
}

export default defineComponent({
  name: "InferenceFormView",
  props: {
    onSubmit: Function as PropType<OnSubmitType>
  },
  computed: {
    pdbFileRule() {
      return this.pdbFile !== null;
    }
  },
  data() : DataType {
    return {
      contig: '',
      hotspots: '',
      numberOfDesigns: 2,
      timesteps: 50,
      pdbFile: null
    }
  },
  methods: {
    async _onSubmit(){
      await this.onSubmit!(this.contig, this.hotspots, this.numberOfDesigns, this.timesteps, this.pdbFile!);
    }
  }
})
</script>

<template>
  <q-form @submit="_onSubmit" class="q-gutter-md">
    <q-input filled v-model="contig" label="Contig" lazy-rules
             :rules="[val => val && val.length > 0 || 'Please type something']">
      <q-tooltip class="text-body1" :offset="[10, 10]" max-width="500px">
        Confirm your chains have the residue numbers you're looking to diffuse over. 100 - Diffuses
        a monomer 100 residues long. 50-100 - Diffuses a hetero-oligomer of lengths 50 and 100.
        5-15/A10-25/30-40 - Builds 5-15 residues N-terminally of A10-25 from the input pdb, followed
        by 30-40 residues to its C-terminus. B1-100/0 100-100 - Generates 100 residue long binders
        to residues 1-100 of chain B.
      </q-tooltip>
    </q-input>
    <q-input filled v-model="hotspots" label="Hotspots">
      <q-tooltip class="text-body1" :offset="[10, 10]" max-width="500px">
        The model optionally readily learns that it should be making an interface which involving
        these hotspot residues. Input is ChainResidueNumber: A100 for residue 100 on chain A.
      </q-tooltip>
    </q-input>
    <q-input filled type="number" v-model="numberOfDesigns" label="Number of desings"
             lazy-rules :rules="[val => val && val > 0 || 'Please type something']">
      <q-tooltip class="text-body1" :offset="[10, 10]" max-width="500px">
        Number of designs to generate
      </q-tooltip>
    </q-input>
    <q-input filled type="number" v-model="timesteps" label="Timesteps" lazy-rules
             :rules="[val => val && val > 0 || 'Please type something']">
      <q-tooltip class="text-body1" :offset="[10, 10]" max-width="500px">
        Desired iterations to generate structure.
      </q-tooltip>
    </q-input>
    <q-file filled bottom-slots accept=".pdb" :rules="[pdbFileRule]" v-model="pdbFile"
            label=".pdb file" counter>
      <template v-slot:prepend>
        <q-icon name="cloud_upload" @click.stop.prevent/>
      </template>
      <template v-slot:append>
        <q-icon name="close" @click.stop.prevent="pdbFile = null" class="cursor-pointer"/>
      </template>

      <template v-slot:hint>
        Upload .pdb file with protein structure
      </template>
    </q-file>
    <div>
      <q-btn label="Run inference" size="large" type="submit" color="positive"/>
    </div>
  </q-form>
</template>
