<script lang="ts">
import {defineComponent, PropType} from 'vue'
import {AminoAcidJobProperties} from "src/features/aminoAcid/types";

interface OnSubmitType {
  (data: AminoAcidJobProperties): Promise<void>;
}

export default defineComponent({
  name: "AminoAcidInferenceForm",
  props: {
    onSubmit: {
      type: Function as PropType<OnSubmitType>,
      required: true
    },
    properties: {
      type: Object as PropType<AminoAcidJobProperties>,
      required: true
    },
    multiple: {
      type: Boolean,
      required: false,
      default: true
    }
  },
  computed: {
    inputLabel() {
      return this.multiple ? ".fasta files (multiple)" : "fasta file"
    },
    hintLabel(){
      return this.multiple ? "Upload .fasta files with amino acid sequences" : "Upload .fasta file with amino acid sequences"
    }
  },
  data(): AminoAcidJobProperties {
    return {
      fastas: this.properties.fastas
    }
  },
  methods: {
    async _onSubmit(evt: Event) {
      await this.onSubmit!(this.$data);
    },
    validateAminoAcid(){
      return this.$data.fastas.length > 0;
    }
  },
})
</script>

<template>
  <q-form @submit="_onSubmit" class="q-gutter-md">
    <q-file filled :multiple="multiple" bottom-slots accept=".fasta" v-model="fastas"
            :rules="[validateAminoAcid]"
            :label="inputLabel" counter>
      <template v-slot:prepend>
        <q-icon name="cloud_upload" @click.stop.prevent/>
      </template>
      <template v-slot:append>
        <q-icon name="close" class="cursor-pointer"/>
      </template>

      <template v-slot:hint>
        {{ hintLabel }}
      </template>
    </q-file>
    <div>
      <q-btn label="Run inference" size="large" type="submit" color="info"/>
    </div>
  </q-form>
</template>
