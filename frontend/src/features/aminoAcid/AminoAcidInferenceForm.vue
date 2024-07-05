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
      <q-btn label="Run inference" size="large" type="submit" color="info"/>
    </div>
  </q-form>
</template>
