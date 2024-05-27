<template>
  <q-card>
    <node-content-header
      name="EsmFold"
      description="Predicts 3D structures from protein sequences"
      :inputs="['FASTA']"
      :outputs="['PDB']"
      imageSrc="/folding_pic.png"
    />
    <node-content-body :results="results" />
  </q-card>
</template>

<script lang="ts">
import { defineComponent } from 'vue';
import NodeContentHeader from 'src/features/drug_discovery/components/workflow/nodeTemplates/NodeContentHeader.vue';
import NodeContentBody from 'src/features/drug_discovery/components/workflow/nodeTemplates/NodeContentBody.vue';

export default defineComponent({
  name: 'EsmFoldNodeContent',
  components: {
    NodeContentHeader,
    NodeContentBody,
  },
  props: {
    results: Array,
  },
});
</script>
