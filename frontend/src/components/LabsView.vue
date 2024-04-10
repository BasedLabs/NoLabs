<template>
  <div class="q-pa-md">
    <!-- Logo at the top, centered -->
    <div class="q-pb-lg flex justify-center">
      <img src="/logo.png" alt="Logo" style="width: 20vw;">
    </div>

    <!-- Vertical list of cards -->
    <div class="q-gutter-md">
      <div v-for="lab in labs" :key="lab.name" class="q-mb-md">
        <q-card class="flex q-pa-md row" flat bordered>
          <q-img class="q-mr-md" width="102px" height="102px" :src="lab.imagePath" alt="Lab Picture"/>
          <div class="col">
            <div class="text-h6 q-pt-md q-mb-sm">{{ lab.name }}</div>
            <div class="q-mb-md">{{ lab.description }}</div>
          </div>
          <div class="column items-center">
            <div class="col">
              <q-btn class="q-ma-lg q-pa-md" :label="'Launch ' + lab.name" color="info"
                     @click="() => launchLab(lab.routeName)"/>
            </div>
          </div>
        </q-card>
      </div>
    </div>
  </div>
</template>

<script lang="ts">
import {defineComponent} from 'vue'

export default defineComponent({
  name: "LabsView",
  data() {
    return {
      labs: [
        {
          name: 'Drug Discovery Lab',
          description: 'Fold proteins, predict binding pockets and run docking',
          routeName: 'Drug discovery experiments',
          imagePath: '/docking_pic.png'
        },
        {
          name: 'Folding Lab',
          description: 'Fold proteins with the latest state-of-the-art models',
          routeName: 'Folding experiments',
          imagePath: '/folding_pic.png'
        },
        {
          name: 'Gene Ontology Lab',
          description: 'Predict gene ontology of proteins',
          routeName: 'Gene ontology experiments',
          imagePath: '/gene_ontology_pic.png'
        },
        {
          name: 'Localisation Lab',
          description: 'Predict localisation of proteins',
          routeName: 'Localisation experiments',
          imagePath: '/localisation_pic.png'
        },
        {
          name: 'Solubility Lab',
          description: 'Predict solubility of proteins',
          routeName: 'Solubility experiments',
          imagePath: '/solubility_pic.png'
        },
        {
          name: 'Conformations Lab',
          description: 'Run simulations using GROMACS',
          routeName: 'Conformations experiments',
          imagePath: '/conformations_pic.png'
        },
        {
          name: 'Protein Design Lab',
          description: 'Design binders using RF-Diffusion',
          routeName: 'Protein design experiments',
          imagePath: '/rf_diffusion_pic.png'
        },
        {
          name: 'Small Molecules Design Lab',
          description: 'Design small molecules using a protein target',
          routeName: 'Small molecules design experiments',
          imagePath: '/small_molecule.jpg'
        }
      ]
    }
  },
  methods: {
    launchLab(routeName: string) {
      this.$router.push({name: routeName});
    }
  }
})
</script>


<style scoped>

</style>
