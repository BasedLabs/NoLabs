<script>
import store from '../storage';
import { api } from '../storage';

export default {
}
</script>

<template>
    <div class="container-fluid row justify-content-center">
      <!-- Amino Acid Lab -->
      <div class="col-md-2 text-center"></div>
      <div class="col-md-2 text-center">
        <RouterLink type="button" class="btn btn-primary btn-md" to="asd">Amino acid lab</RouterLink>
      </div>

      <!-- Drug Discovery Lab -->
      <div class="col-md-2 text-center">
        <RouterLink type="button" class="btn btn-primary btn-md" to="/drug-target-lab">Drug discovery lab</RouterLink>
      </div>

      <!-- Conformations lab Lab -->
      <div class="col-md-2 text-center">
        <RouterLink type="button" class="btn btn-primary btn-md" to="/conformations">Conformations lab (pre-alpha)</RouterLink>
      </div>

      <!-- Conformations lab Lab -->
      <div class="col-md-2 text-center">
        <RouterLink type="button" class="btn btn-primary btn-md" to="/protein-design">Protein design lab</RouterLink>
      </div>
      <div class="col-md-2 text-center"></div>
    </div>
</template>