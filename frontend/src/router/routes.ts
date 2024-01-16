const routes = [
  {
    path: '/labs',
    component: () => import('../layouts/MainLayout.vue'),
    children: [
      {
        path: 'solubility',
        component: () => import('../features/solubility/Page.vue')
      },
      {
        path: 'drug-discovery',
        component: () => import('../features/drug_discovery/Page.vue')
      },
      {
        path: 'drug-discovery/experiment/:experimentId',
        component: () => import('../features/drug_discovery/ExperimentViewer.vue'),
        name: 'ExperimentViewer', 
        props: true 
      }
    ]
  }
];

export default routes
