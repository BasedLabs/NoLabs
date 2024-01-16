<template>
  <div class="q-pa-md">
    <q-layout view="lHh Lpr lFf">
      <q-header class="bg-black">
        <q-toolbar>
          <q-toolbar-title>
            !Labs
          </q-toolbar-title>

          <q-btn flat dense round icon="home" aria-label="Home" />
        </q-toolbar>
      </q-header>
      <q-page-container>
        <router-view />
      </q-page-container>
    </q-layout>
  </div>
</template>

<script lang="ts">
import { defineComponent } from 'vue';

export default defineComponent({
  name: 'MainLayout',
  data() {
    return {
      drawer: false
    }
  },
  components: {

  }
})
</script>
