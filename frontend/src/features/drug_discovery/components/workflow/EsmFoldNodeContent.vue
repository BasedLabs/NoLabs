<template>
  <q-card>
    <q-card-section class="bg-black row items-center q-mt-sm text-white rounded-borders q-pa-md">
      <div class="col-4">
        <q-img class="rounded-borders" width="102px" height="102px" :src="'/folding_pic.png'" alt="Lab Picture"/>
      </div>
      <div class="col-8 flex flex-column justify-center text-center q-mb-lg">
        <h5 class="q-mb-sm">EsmFold</h5>
        <h8>Predicts 3D structures from protein sequences</h8>
      </div>
    </q-card-section>

    <!-- List of jobs -->
    <q-card-section class="job-section">
      <h7>Jobs</h7>
      <q-list v-if="jobs && jobs.length">
        <q-item v-for="job in jobs" :key="job.id">
          <q-item-section>{{ job.name }}</q-item-section>
        </q-item>
      </q-list>
    </q-card-section>

    <!-- List of results -->
    <q-card-section>
      <h7>Results</h7>
      <q-list v-if="results && results.length">
        <q-item v-for="result in results" :key="result.id">
          <q-item-section>{{ result.name }}</q-item-section>
        </q-item>
      </q-list>
    </q-card-section>
  </q-card>
</template>

<script>
import { defineComponent } from 'vue';

export default defineComponent({
  name: 'EsmFoldNodeContent',
  props: {
    jobs: Array,
    results: Array,
  },
});
</script>

<style scoped>
.job-section {
  border: 1px dashed #ccc;
  padding: 10px;
}
</style>
