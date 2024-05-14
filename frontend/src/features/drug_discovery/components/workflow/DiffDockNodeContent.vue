<template>
  <q-card class="q-pa-md">
    <q-card-section class="bg-black row no-wrap items-center q-mt-sm text-white rounded-borders">
      <q-img class="rounded-borders" width="102px" height="102px" :src="'/docking_pic.png'" alt="Lab Picture"/>
      <h5 class="q-pl-md">DiffDock</h5>
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
      <q-list  v-if="results && results.length">
        <q-item v-for="result in results" :key="result.id">
          <q-item-section>{{ result.name }}</q-item-section>
        </q-item>
      </q-list>
    </q-card-section>
  </q-card>
</template>

<script>
import {defineComponent} from 'vue';

export default defineComponent({
  name: 'DiffDockNodeContent',
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
