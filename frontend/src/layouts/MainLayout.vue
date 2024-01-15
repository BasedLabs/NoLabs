<template>
  <div class="q-pa-md">
    <q-layout view="lHh Lpr lFf">
      <q-header class="bg-black">
        <q-toolbar>
          <q-btn flat @click="drawer = !drawer" round dense icon="menu" />
          <q-toolbar-title>
            !Labs
          </q-toolbar-title>

          <q-btn flat dense round icon="home" aria-label="Home" />
        </q-toolbar>
      </q-header>

      <q-drawer v-model="drawer" show-if-above :width="200" :breakpoint="500">
        <q-scroll-area class="fit">
          <ExperimentsList />
        </q-scroll-area>
      </q-drawer>

      <q-page-container>
        <router-view />
      </q-page-container>
    </q-layout>
  </div>
</template>

<script lang="ts">
import { defineComponent } from 'vue';
import ExperimentsList from 'components/ExperimentsList.vue';

export default defineComponent({
  name: 'MainLayout',
  data() {
    return {
      drawer: false
    }
  },
  components: {
    ExperimentsList
  }
})
</script>
