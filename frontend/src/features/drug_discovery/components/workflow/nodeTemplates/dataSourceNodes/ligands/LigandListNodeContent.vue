<template>
  <q-card>
    <LigandsList v-if="experimentId" :experiment-id="experimentId" />
  </q-card>
</template>

<script lang="ts">
import { defineComponent } from 'vue';
import LigandsList from './LigandsList.vue';

export default defineComponent({
  name: 'LigandsNodeContent',
  components: { LigandsList },
  props: {
    experimentId: String,
  },
  data() {
    return {
      loading: true,
    };
  },
});
</script>