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
        component: () => import('../features/drug_discovery/ExperimentsList.vue')
      },
      {
        path: 'drug-discovery/experiment/:experimentId',
        component: () => import('../features/drug_discovery/ExperimentNavigation.vue'),
        name: 'ExperimentNavigation',
        props: true,
        children: [
          {
            path: 'drug-discovery/experiment/:experimentId/setup',
            component: () => import('../features/drug_discovery/components/ExperimentSetup.vue'),
            name: 'ExperimentSetup',
            props: true
          },
          {
            path: 'drug-discovery/experiment/:experimentId/running-jobs',
            component: () => import('../features/drug_discovery/components/RunningJobs.vue'),
            name: 'RunningJobs',
            props: true
          },
        ]
      }
    ]
  }
];

export default routes
