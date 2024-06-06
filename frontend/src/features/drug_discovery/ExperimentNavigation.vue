<template>
  <q-page class="bg-black q-pl-md">
    <BioBuddyChat v-if="bioBuddyEnabled" :experiment-id="experiment.experimentId as string" />
    <div class="row no-wrap items-center">
      <q-btn-dropdown color="primary" label="Add Node" icon="add" dense persistent>
        <q-list>
          <q-item v-for="option in componentOptions" :key="option.name" clickable v-close-popup
            @click="addComponent(option)">
            <q-item-section>
              <q-item-label>{{ option.name }}</q-item-label>
              <q-item-label caption>{{ option.description }}</q-item-label>
            </q-item-section>
          </q-item>
        </q-list>
      </q-btn-dropdown>
      <q-space />
      <q-btn class="q-ma-sm" color="green" icon="not_started" @click="startWorkflow">Start workflow</q-btn>
    </div>

    <div class="map-container">
      <VueFlow class="workflow" v-if="elements" :nodes="elements.nodes" @nodeDragStop="onNodeDragStopHandler"
        :edges="elements.edges" @connect="onConnect" fit-view-on-init>
        <template #node-Proteins="{ id }">
          <ProteinListNode :nodeId="id" :onDeleteNode="onDeleteNode" :onOpenSettings="openSettings"
            :onOpenDialog="openNodeDialog" />
        </template>
        <template #node-ligand-list="{ id }">
          <LigandListNode :nodeId="id" :onDeleteNode="onDeleteNode" :onOpenSettings="openSettings"
            :onOpenDialog="openNodeDialog" />
        </template>
        <template #node-custom="{ id }">
          <JobNode :nodeId="id" :onDeleteNode="onDeleteNode" :onOpenSettings="openSettings"
            :onOpenDialog="openNodeDialog" />
        </template>
        <template #edge-custom="customEdgeProps">
          <ErrorEdge :id="customEdgeProps.id" :source-x="customEdgeProps.sourceX" :source-y="customEdgeProps.sourceY"
            :target-x="customEdgeProps.targetX" :target-y="customEdgeProps.targetY"
            :source-position="customEdgeProps.sourcePosition" :target-position="customEdgeProps.targetPosition"
            :data="customEdgeProps.data" :marker-end="customEdgeProps.markerEnd" :style="customEdgeProps.style" />
        </template>
      </VueFlow>
    </div>
  </q-page>

  <q-drawer v-model="sideMenuOpen" v-show="sideMenuOpen" bordered content-style="background-color: white;" side="right">
    <!-- Close button -->
    <q-btn @click="closeSideMenu" class="q-ma-md" flat round dense icon="close" />

    <!-- Side menu content -->
    <q-card>
      <q-card-section>
        <q-item>
          <q-item-section>
            <q-item-label v-if="selectedNode">{{ selectedNode.label }}</q-item-label>
            <q-item-label v-if="selectedNode" caption>{{ selectedNode.description }}</q-item-label>
          </q-item-section>
        </q-item>
      </q-card-section>
    </q-card>
  </q-drawer>

  <q-dialog v-model="modalOpen" persistent>
    <q-card style="min-width: 70vw; min-height: 70vh;">
      <q-card-actions align="right">
        <q-btn flat round dense icon="close" v-close-popup @click="closeModal" />
      </q-card-actions>
      <q-card-section>
        <ProteinListNodeContent v-if="experiment.experimentId && selectedNode && selectedNode.type === 'protein-list'"
          :experiment-id="experiment.experimentId" />
        <LigandListNodeContent v-if="experiment.experimentId && selectedNode && selectedNode.type === 'ligand-list'"
          :experiment-id="experiment.experimentId" />
        <JobNodeContent v-if="experiment.experimentId && selectedNode && selectedNode.type === 'esmfold'"
          :experiment-id="experiment.experimentId" :name="selectedNode?.name" :description="selectedNode?.description"
          :jobIds="selectedNode?.data.jobIds" />
      </q-card-section>
    </q-card>
  </q-dialog>
</template>

<script lang="ts">
import '@vue-flow/core/dist/style.css';
import '@vue-flow/core/dist/theme-default.css';

import BioBuddyChat from "src/features/biobuddy/BioBuddyChat.vue";
import ProteinListNode from "./components/workflow/nodeTemplates/dataSourceNodes/ProteinListNode.vue";
import ProteinListNodeContent from "./components/workflow/nodeTemplates/dataSourceNodes/ProteinListNodeContent.vue";
import LigandListNode from "./components/workflow/LigandListNode.vue";
import LigandListNodeContent from "./components/workflow/LigandListNodeContent.vue";
import ErrorEdge from "./components/workflow/nodeTemplates/ErrorEdge.vue"
import { useDrugDiscoveryStore } from "./storage";
import { defineComponent } from "vue";
import { Edge, Node as FlowNode } from '@vue-flow/core';
import { VueFlow } from '@vue-flow/core';
import JobNode from "./components/workflow/JobNode.vue";
import JobNodeContent from "./components/workflow/JobNodeContent.vue";
import { startWorkflowforExperiment, checkBiobuddyEnabled } from 'src/features/drug_discovery/refinedApi';
import { useWorkflowStore } from 'src/features/drug_discovery/components/workflow/storage';

// Define custom Node type
interface Node extends FlowNode {
  id: string;
  name: string;
  type: string;
  inputs: string[];
  outputs: string[];
  jobIds: string[];
  description: string;
  error: string;
}

export default defineComponent({
  name: "ExperimentNavigation",
  components: {
    BioBuddyChat,
    VueFlow,
    ProteinListNode,
    ProteinListNodeContent,
    LigandListNode,
    LigandListNodeContent,
    JobNode,
    JobNodeContent,
    ErrorEdge
  },
  data() {
    return {
      step: 1, // This can be a reactive property based on the current route if needed
      experiment: {
        experimentId: null as string | null,
        metadata: null
      },
      elements: {
        nodes: [] as Node[],
        edges: [] as Edge[]
      },
      allowedTypes: ["Proteins", "Ligands", "DNA"],
      sideMenuOpen: false,
      modalOpen: false,
      selectedNode: null as Node | null,
      selectedComponent: null as string | null,
      specialNodeProps: [],
      splitterModel: 20,
      bioBuddyEnabled: false,
      workflowId: "19fd2fa8-9a01-4755-adb2-b11a6b302e7f" // Example workflow ID
    };
  },
  computed: {
    componentOptions() {
      const workflowStore = useWorkflowStore();
      return workflowStore.componentOptions;
    }
  },
  async mounted() {
    const store = useDrugDiscoveryStore();
    this.experiment.experimentId = this.$route.params.experimentId as string;
    this.experiment.metadata = await store.getExperimentMetaData(this.experiment.experimentId);
    try {
      this.bioBuddyEnabled = await checkBiobuddyEnabled();
    } catch (error) {
      console.error('Error checking BioBuddy enabled status:', error);
      this.bioBuddyEnabled = false;
    }

    const workflowStore = useWorkflowStore();
    await workflowStore.fetchWorkflow(this.workflowId);
    this.elements = workflowStore.elements;

    setInterval(() => {
      workflowStore.pollWorkflow();
    }, 2000); // Poll every 2 seconds
  },
  methods: {
    addComponent(option: any) {
      const workflowStore = useWorkflowStore();
      workflowStore.addComponent(option);
    },
    onConnect(params: { source: string, target: string, sourceHandle: string, targetHandle: string }) {
      const workflowStore = useWorkflowStore();
      workflowStore.onConnect(params);
    },
    onDeleteNode(nodeId: string) {
      const workflowStore = useWorkflowStore();
      workflowStore.onDeleteNode(nodeId);
    },
    openSettings(nodeId: string) {
      const workflowStore = useWorkflowStore();
      const node = workflowStore.getNodeById(nodeId);
      if (node) {
        this.selectedNode = node;
        this.sideMenuOpen = true;
      }
    },
    closeSideMenu() {
      this.sideMenuOpen = false;
    },
    closeModal() {
      this.modalOpen = false;
    },
    startWorkflow() {
      startWorkflowforExperiment(this.experiment.experimentId as string);
    },
    onNodeDragStopHandler(event: { node: Node }) {
      const workflowStore = useWorkflowStore();
      workflowStore.onNodeDragStopHandler(event);
    },
    async onExperimentNameChange(newExperimentName: string) {
      const store = useDrugDiscoveryStore();
      await store.changeExperimentName(this.experiment.experimentId as string, newExperimentName);
    },
    openNodeDialog(nodeId: string) {
      const workflowStore = useWorkflowStore();
      const node = workflowStore.getNodeById(nodeId);
      if (node) {
        this.selectedNode = node;
        this.modalOpen = true;
      }
    }
  }
})
</script>

<style scoped>
body {
  overflow-x: hidden;
}

.map-container {
  width: 100vw;
  height: 100vh;
}

.workflow {
  height: 100%;
  width: 100%;
}

.vue-flow__edge {
  z-index: 10 !important;
  position: relative;
}

.vue-flow__node {
  z-index: 5 !important;
  position: relative;
}
</style>
