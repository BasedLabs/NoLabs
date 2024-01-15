const routes = [
  {
    path: '/labs',
    component: () => import('layouts/MainLayout.vue'),
    children: [
      { path: 'solubility', component: () => import('../features/solubility/Page.vue') }
    ]
  }
]

export default routes
