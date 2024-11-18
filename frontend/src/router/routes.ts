const experimentsRoutes = [
  {
    path: 'experiments',
    component: () => import('../features/workflow/ExperimentsView.vue'),
    name: 'Experiments'
  },
]

const workflowRoutes = [
  {
    path: 'workflow/:experimentId',
    component: () => import('../features/workflow/WorkflowView.vue'),
    name: 'Workflow view'
  },
]

const proteinDesignRoutes = [
  {
    path: 'rfdiffusion/job/:jobId',
    component: () => import('../features/proteinDesign/ProteinDesignJobView.vue'),
    name: 'Rfdiffusion',
    props: true
  },
];

const foldingRoutes = [
  {
    path: 'folding/job/:jobId',
    component: () => import('src/features/aminoAcid/folding/FoldingJobView.vue'),
    name: 'Folding',
    props: true
  },
];

const proteinMpnnRoutes = [
  {
    path: 'proteinMpnn/job/:jobId',
    component: () => import('../features/workflow/components/jobs/ProteinMpnnJob.vue'),
    name: 'ProteinMPNN',
    props: true
  }
]

const routes = [
  {
    path: '/',
    name: 'Main',
    component: () => import('src/layouts/MainLayout.vue'),
    children: [
      ...experimentsRoutes,
      ...workflowRoutes,
      ...proteinDesignRoutes,
      ...foldingRoutes,
      ...proteinMpnnRoutes
    ]
  },
  {
    path: '/:catchAll(.*)*',
    component: () => import('pages/ErrorNotFound.vue'),
  },
];

export default routes
