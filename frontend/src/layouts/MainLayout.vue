<template>
  <div>
    <q-layout view="lHh Lpr lFf">
      <q-header class="bg-black">
        <q-toolbar>
          <q-toolbar-title>
            <router-link :to="{ name: 'Experiments' }">
              <q-img src="/logo.png" fit="scale-down" style="max-width: 100px">
              </q-img>
            </router-link>
          </q-toolbar-title>
        </q-toolbar>
      </q-header>
      <q-page-container>
        <router-view />
      </q-page-container>
    </q-layout>
  </div>
</template>

<script lang="ts">
import { defineComponent } from "vue";

export default defineComponent({
  name: "MainLayout",
  mounted(){
    debugger;
    if(this.$route.name === 'Main'){
      this.$router.push({name: 'Experiments'});
    }
  }
});
</script>
