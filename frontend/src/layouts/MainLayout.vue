<template>
  <div>
    <q-layout view="lHh Lpr lFf">
      <q-header class="bg-black">
        <q-toolbar>
          <q-toolbar-title>
            <q-img src="/logo.png" fit="scale-down" style="max-width: 100px;">
            </q-img>
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
