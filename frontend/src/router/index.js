import { createRouter, createWebHistory } from 'vue-router'
import AminoAcidLabView from '../views/AminoAcidLab.vue';
import DrugTargetLabView from '../views/DrugTargetLab.vue';
import ConformationsLab from '../views/ConformationsLab.vue';
import ProteinDesignLab from '../views/ProteinDesignLab.vue';
import ProteinViewerView from '../views/ProteinViewerView.vue';

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
    {
      path: '/conformations',
      name: 'conformationsLab',
      component: ConformationsLab
    }
    ,
    {
      path: '/protein-design',
      name: 'proteinDesign',
      component: ProteinDesignLab
    },
    {
      path: '/protein-viewer',
      name: 'proteinViewer',
      component: ProteinViewerView
    }
  ]
})

export default router
