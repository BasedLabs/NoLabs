import { createRouter, createWebHistory } from 'vue-router'
import AminoAcidLabView from '../views/AminoAcidLab.vue';
import DrugTargetLabView from '../views/DrugTargetLab.vue';

const router = createRouter({
  history: createWebHistory(import.meta.env.BASE_URL),
  routes: [
    {
      path: '/amino-acid-lab',
      name: 'aminoAcid',
      component: AminoAcidLabView
    },
    {
      path: '/drug-target-lab',
      name: 'drugTargetLab',
      component: DrugTargetLabView
    },
  ]
})

export default router
