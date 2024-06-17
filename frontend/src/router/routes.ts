const drugDiscoveryRoutes = [
  {
    path: 'drug-discovery',
    name: 'Drug discovery experiments',
    component: () => import('src/features/drug_discovery/ExperimentsView.vue')
  },
  {
    path: 'drug-discovery/experiment/:experimentId',
    component: () => import('src/features/drug_discovery/ExperimentNavigation.vue'),
    name: 'Drug discovery',
    props: true,
    children: [
      {
        path: '/drug-discovery/experiment/:experimentId/upload-targets',
        component: () => import('../features/drug_discovery/components/UploadTargetsStep.vue'),
        name: 'Upload targets',
        props: true
      },
      {
        path: '/drug-discovery/experiment/:experimentId/upload-ligands',
        component: () => import('../features/drug_discovery/components/UploadLigandsStep.vue'),
        name: 'Upload ligands',
        props: true
      },
      {
        path: '/drug-discovery/experiment/:experimentId/run-docking',
        component: () => import('../features/drug_discovery/components/RunDockingStep.vue'),
        name: 'Run docking',
        props: true
      },
    ]
  }
];


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

const labsRoute = [
  {
    path: 'labs',
    name: 'Labs',
    component: () => import("src/components/LabsView.vue")
  },
];

const routes = [
  {
    path: '/',
    name: 'Main',
    component: () => import('src/layouts/MainLayout.vue'),
    children: [
      ...labsRoute,
      ...drugDiscoveryRoutes,
      ...proteinDesignRoutes,
      ...conformationsRoutes,
      ...localisationRoutes,
      ...solubilityRoutes,
      ...geneOntologyRoutes,
      ...foldingRoutes,
      ...smallMoleculeDesignRoutes
    ]
  }
];

export default routes
