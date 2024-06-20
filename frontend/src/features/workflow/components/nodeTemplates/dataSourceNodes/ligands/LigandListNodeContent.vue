<template>
  <q-card>
    <LigandsList v-if="experimentId" :experiment-id="experimentId" :nodeId="nodeId" />
  </q-card>
</template>

<script lang="ts">
import { defineComponent } from 'vue';
import LigandsList from './LigandsList.vue';

export default defineComponent({
  name: 'LigandsNodeContent',
  components: { LigandsList },
  props: {
    experimentId: {
      type: String,
      required: true,
    },
    nodeId: {
      type: String,
      required: true,
    },
  },
  data() {
    return {
      loading: true,
    };
  },
});
</script>