<template>
  <q-card>
    <NodeContentHeader :name="name" :description="description" :inputs="['FASTA']" :outputs="['PDB']"
      imageSrc="/folding_pic.png" />
    <NodeContentBody :nodeId="nodeId" :name="name" />
  </q-card>
</template>

<script lang="ts">
import { defineComponent } from 'vue';
import NodeContentHeader from 'src/features/drug_discovery/components/workflow/nodeTemplates/NodeContentHeader.vue';
import NodeContentBody from 'src/features/drug_discovery/components/workflow/nodeTemplates/NodeContentBody.vue';

export default defineComponent({
  name: 'JobNodeContent',
  components: {
    NodeContentHeader,
    NodeContentBody,
  },
  props: {
    nodeId: {
      type: String,
      required: true
    },
    name: {
      type: String,
      required: true
    },
    description: {
      type: String,
      required: true
    }
  },
  mounted() {
    console.log(this.name)
  }
});
</script>
