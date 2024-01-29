<script lang="ts">
import {defineComponent} from 'vue'

export default defineComponent({
  name: "LabsView"
})
</script>

<template>
  <div class="q-pa-md row items-start q-gutter-md justify-center">
    <q-card flat bordered class="my-card">
      <q-card-section>
        <q-btn size="xl" color="positive"  label="Folding lab"
               @click="() => {this.$router.push({name: 'Folding'})}"/>
      </q-card-section>
    </q-card>
    <q-card flat bordered class="my-card">
      <q-card-section>
        <q-btn size="xl" color="positive"  label="Gene ontology lab"
               @click="() => {this.$router.push({name: 'GeneOntology'})}"/>
      </q-card-section>
    </q-card>
    <q-card flat bordered class="my-card">
      <q-card-section>
        <q-btn size="xl" color="positive"  label="Localisation lab"
               @click="() => {this.$router.push({name: 'Localisation'})}"/>
      </q-card-section>
    </q-card>
    <q-card flat bordered class="my-card">
      <q-card-section>
        <q-btn size="xl" color="positive"  label="Solubility lab"
               @click="() => {this.$router.push({name: 'Solubility'})}"/>
      </q-card-section>
    </q-card>
    <q-card flat bordered class="my-card">
      <q-card-section>
        <q-btn size="xl" color="positive"  label="Conformations lab"
               @click="() => {this.$router.push({name: 'Conformations'})}"/>
      </q-card-section>
    </q-card>
    <q-card flat bordered class="my-card">
      <q-card-section>
        <q-btn size="xl" color="positive"  label="Drug discovery lab"
               @click="() => {this.$router.push({name: 'DrugDiscovery'})}"/>
      </q-card-section>
    </q-card>
    <q-card flat bordered class="my-card">
      <q-card-section>
        <q-btn size="xl" color="positive"  label="Protein design lab"
               @click="() => {this.$router.push({name: 'ProteinDesign'})}"/>
      </q-card-section>
    </q-card>
  </div>
</template>

<style scoped>

</style>