<script setup>
import { RouterLink, RouterView } from 'vue-router';
</script>

<template>
  <div class="text-center" v-if="this.$route.path === '/'">
    <RouterLink type="button" class="btn btn-primary btn-md" to="/amino-acid-lab">Amino acid lab</RouterLink>
    <RouterLink type="button" class="btn btn-primary btn-md ms-1" to="/drug-target-lab">Drug target discovery lab
    </RouterLink>
  </div>
  <RouterView />
</template>