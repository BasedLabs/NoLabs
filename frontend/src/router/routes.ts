const solubilityRoutes = [
  {
    path: 'solubility',
    component: () => import('../features/solubility/Page.vue')
  },
]

const drugDiscoveryRoutes = [
  {
    path: 'drug-discovery',
    component: () => import('../features/drug_discovery/Page.vue')
  },
  {
    path: 'drug-discovery/experiment/:experimentId',
    component: () => import('../features/drug_discovery/ExperimentViewer.vue'),
    name: 'ExperimentViewer',
    props: true
  },
];

const proteinDesignRoutes = [
  {
    path: 'protein-design',
    component: () => import('../features/proteinDesign/ExperimentsView.vue')
  },
  {
    path: 'protein-design/experiment/:experimentId',
    component: () => import('../features/proteinDesign/ExperimentView.vue'),
    name: 'ProteinDesignExperimentView',
    props: true
  },
];

const routes = [
  {
    path: '/labs',
    component: () => import('../layouts/MainLayout.vue'),
    children: [
      ...solubilityRoutes,
      ...drugDiscoveryRoutes,
      ...proteinDesignRoutes
    ]
  }
];

export default routes
