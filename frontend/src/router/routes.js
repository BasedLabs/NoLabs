const drugDiscoveryRoutes = [
    {
        path: 'drug-discovery',
        component: () => import('src/features/drug_discovery/Page.vue')
    },
    {
        path: 'drug-discovery/experiment/:experimentId',
        component: () => import('src/features/drug_discovery/ExperimentViewer.vue'),
        name: 'ExperimentViewer',
        props: true
    },
];


const proteinDesignRoutes = [
    {
        path: 'protein-design',
        component: () => import("src/features/proteinDesign/ExperimentsView.vue")
    },
    {
        path: 'protein-design/experiment/:experimentId',
        component: () => import('src/features/proteinDesign/ExperimentView.vue'),
        name: 'ProteinDesignExperimentView',
        props: true
    },
];

const conformationsRoutes = [
    {
        path: 'conformations',
        component: () => import("src/features/conformations/ExperimentsView.vue")
    },
    {
        path: 'conformations/experiment/:experimentId',
        component: () => import('src/features/conformations/ExperimentView.vue'),
        name: 'ConformationsExperimentView',
        props: true
    }
]

const routes = [
    {
        path: '/labs',
        name: 'Labs',
        component: () => import('src/layouts/MainLayout.vue'),
        children: [
            ...drugDiscoveryRoutes,
            ...proteinDesignRoutes,
            ...conformationsRoutes
        ]
    }
];

export default routes
