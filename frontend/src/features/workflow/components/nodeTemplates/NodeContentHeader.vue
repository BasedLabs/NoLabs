<template>
  <q-card-section class="bg-black row items-center q-mt-sm text-white rounded-borders q-pl-md q-pr-md">
    <div class="col-4 text-right">
      <q-img class="rounded-borders" width="100px" height="100px" :src="imageSrc" alt="Component Picture"/>
    </div>
    <div class="col-8 flex column justify-center text-center q-mb-md">
      <h5 class="q-mb-sm">{{ name }}</h5>
      <div class="q-pa-sm">{{ description }}</div>
    </div>
  </q-card-section>
</template>

<script lang="ts">
import {defineComponent} from 'vue';

export default defineComponent({
  name: 'NodeContentHeader',
  props: {
    name: {
      type: String,
      required: true
    },
    description: {
      type: String,
      required: true
    },
    inputs: {
      type: Array,
      required: true
    },
    outputs: {
      type: Array,
      required: true
    },
    imageSrc: {
      type: String,
      required: false
    }
  }
});
</script>
