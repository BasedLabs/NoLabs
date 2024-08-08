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
    path: 'protein-design/job/:jobId',
    component: () => import('../features/proteinDesign/ProteinDesignJobView.vue'),
    name: 'Protein design',
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

const blastRoutes = [
  {
    path: ':experimentId/blast/job/:jobId',
    component: () => import('src/features/workflow/components/jobs/BlastJob.vue'),
    name: 'Blast',
    props: true
  },
];

const localisationRoutes = [
  {
    path: 'localisation/job/:jobId',
    component: () => import('../features/aminoAcid/localisation/LocalisationJobView.vue'),
    name: 'Localisation',
    props: true
  },
];

const geneOntologyRoutes = [
  {
    path: 'gene-ontology/job/:jobId',
    component: () => import('src/features/aminoAcid/geneOntology/GeneOntologyJobView.vue'),
    name: 'Gene ontology',
    props: true
  },
];

const solubilityRoutes = [
  {
    path: 'solubility/job/:jobId',
    component: () => import('../features/aminoAcid/solubility/SolubilityJobView.vue'),
    name: 'Solubility',
    props: true
  },
];

const conformationsRoutes = [
  {
    path: 'conformations/job/:jobId',
    component: () => import('../features/conformations/ConformationsJobView.vue'),
    name: 'Conformations',
    props: true
  }
]

const smallMoleculeDesignRoutes = [
  {
    path: 'small-molecules/job/:jobId',
    component: () => import('../features/smallMoleculeDesign/SmallMoleculesDesignJobView.vue'),
    name: 'Small molecules design',
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
      ...conformationsRoutes,
      ...localisationRoutes,
      ...solubilityRoutes,
      ...geneOntologyRoutes,
      ...foldingRoutes,
      ...blastRoutes,
      ...smallMoleculeDesignRoutes
    ]
  }
];

export default routes
