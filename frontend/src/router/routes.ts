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

const blastRoutes = [
  {
    path: ':experimentId/blast/job/:jobId',
    component: () => import('src/features/workflow/components/jobs/BlastJob.vue'),
    name: 'Blast',
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

const proteinAffinityCharacterizationRoutes = [
  {
    path: 'proteinAffinityCharacterization/job/:jobId',
    component: () => import('../features/workflow/components/jobs/ProteinAffinityCharacterizationJob.vue'),
    name: 'ProteinAffinityCharacterization',
    props: true
  }
]

const proteinBindingScreeningRoutes = [
  {
    path: 'proteinBindingScreening/job/:jobId',
    component: () => import('../features/workflow/components/jobs/ProteinBindingScreeningJob.vue'),
    name: 'ProteinBindingScreening',
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
      ...blastRoutes,
      ...proteinMpnnRoutes,
      ...proteinAffinityCharacterizationRoutes,
      ...proteinBindingScreeningRoutes
    ]
  },
  {
    path: '/:catchAll(.*)*',
    component: () => import('pages/ErrorNotFound.vue'),
  },
];

export default routes
