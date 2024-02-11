<script lang="ts">
import {defineComponent, PropType} from 'vue'
import {AminoAcidExperimentProperties} from "src/features/aminoAcid/types";

interface OnSubmitType {
  (data: AminoAcidExperimentProperties): Promise<void>;
}

export default defineComponent({
  name: "AminoAcidInferenceForm",
  props: {
    onSubmit: {
      type: Function as PropType<OnSubmitType>,
      required: true
    },
    properties: {
      type: Object as PropType<AminoAcidExperimentProperties>,
      required: true
    }
  },
  data(): AminoAcidExperimentProperties {
    return {
      aminoAcidSequence: this.properties.aminoAcidSequence,
      fastas: this.properties.fastas
    }
  },
  methods: {
    async _onSubmit(evt: Event) {
      await this.onSubmit!(this.$data);
    },
    validateAminoAcid(){
      if(this.$data.aminoAcidSequence){
        return true;
      }

      return this.$data.fastas.length > 0;
    }
  },
})
</script>

<template>
  <q-form @submit="_onSubmit" class="q-gutter-md">
    <q-input filled v-model="aminoAcidSequence" label="Amino acid sequence"
    :rules="[validateAminoAcid]">
      <q-tooltip class="text-body1" :offset="[10, 10]" max-width="500px">
        Enter amino acid sequence or fasta
      </q-tooltip>
    </q-input>
    <q-file filled multiple bottom-slots accept=".fasta" v-model="fastas"
            :rules="[validateAminoAcid]"
            label=".fasta files (multiple)" counter>
      <template v-slot:prepend>
        <q-icon name="cloud_upload" @click.stop.prevent/>
      </template>
      <template v-slot:append>
        <q-icon name="close" class="cursor-pointer"/>
      </template>

      <template v-slot:hint>
        Upload .fasta files with amino acid sequences
      </template>
    </q-file>
    <div>
      <q-btn label="Run inference" size="large" type="submit" color="positive"/>
    </div>
  </q-form>
</template>
