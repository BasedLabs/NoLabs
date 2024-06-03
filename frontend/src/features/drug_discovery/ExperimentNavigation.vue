<template>
  <q-page class="bg-black q-pl-md">
    <BioBuddyChat v-if="bioBuddyEnabled" :experiment-id="experiment.experimentId as string" />
    <div class="row no-wrap items-center">
      <button @click="generateWorkflow">Generate Workflow</button>
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
      <q-btn class="q-ma-sm" color="green" icon="not_started">Run workflow</q-btn>
    </div>

    <div class="map-container">
      <VueFlow class="workflow" v-if="elements" :nodes="elements.nodes" :nodeDragStop="onNodeDragStopHandler"
        :edges="elements.edges" @connect="onConnect" fit-view-on-init>
        <template #node-Proteins="{ id }">
          <ProteinListNode :nodeId="id" :inputs="elements.nodes.find(n => n.id === id)?.data.inputs"
            :outputs="elements.nodes.find(n => n.id === id)?.data.outputs" :onDeleteNode="onDeleteNode"
            :onOpenSettings="openSettings" :onOpenDialog="openNodeDialog" />
        </template>
        <template #node-ligand-list="{ id }">
          <LigandListNode :nodeId="id" :inputs="elements.nodes.find(n => n.id === id)?.data.inputs"
            :outputs="elements.nodes.find(n => n.id === id)?.data.outputs" :onDeleteNode="onDeleteNode"
            :onOpenSettings="openSettings" :onOpenDialog="openNodeDialog" />
        </template>
        <template #node-Folding="{ id }">
          <JobNode :nodeId="id" :inputs="elements.nodes.find(n => n.id === id)?.data.inputs"
            :outputs="elements.nodes.find(n => n.id === id)?.data.outputs" 
            :jobIds="elements.nodes.find(n => n.id === id)?.data.jobIds"
            :onDeleteNode="onDeleteNode"
            :onOpenSettings="openSettings" :onOpenDialog="openNodeDialog" />
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
          :experiment-id="experiment.experimentId" :jobIds="selectedNode?.data.jobIds" />
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
import { useDrugDiscoveryStore } from "./storage";
import { defineComponent } from "vue";
import { checkBioBuddyEnabled } from "../biobuddy/api";
import { Edge, Node as FlowNode } from '@vue-flow/core';
import { VueFlow } from '@vue-flow/core';
import JobNode from "./components/workflow/JobNode.vue";
import JobNodeContent from "./components/workflow/JobNodeContent.vue";
import { getWorkflow, sendWorkflowUpdate } from 'src/features/drug_discovery/refinedApi';
import {
  ComponentModel_Output,
  WorkflowSchemaModel_Input,
  ComponentModel_Input,
  PropertyModel_Output,
  WorkflowComponentModel,
  MappingModel
} from 'src/refinedApi/client';

import { v4 as uuidv4 } from 'uuid';

// Define custom Node type
interface Node extends FlowNode {
  id: string;
  name: string;
  type: string;
  inputs: string[];
  outputs: string[];
  jobIds: string[];
  description: string;
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
    JobNodeContent
  },
  data() {
    return {
      step: 1, // This can be a reactive property based on the current route if needed
      experiment: {
        experimentId: null as string | null,
        metadata: null
      },
      sideMenuOpen: false,
      modalOpen: false,
      selectedNode: null as Node | null,
      selectedComponent: null as string | null,
      specialNodeProps: [],
      splitterModel: 20,
      bioBuddyEnabled: false,
      componentOptions: [] as Array<{ name: string; type: string; inputs: Record<string, PropertyModel_Output>, outputs: Record<string, PropertyModel_Output>, description: string }>,
      elements: {
        nodes: [] as Node[],
        edges: [] as Edge[]
      },
      workflowId: "d35dd2d8-8962-4ce2-b0c0-c9cbb5cfee4e" // Example workflow ID
    };
  },
  async mounted() {
    const store = useDrugDiscoveryStore();
    this.experiment.experimentId = this.$route.params.experimentId as string;
    this.experiment.metadata = await store.getExperimentMetaData(this.experiment.experimentId);
    try {
      const response = await checkBioBuddyEnabled();
      this.bioBuddyEnabled = response.enabled;
    } catch (error) {
      console.error('Error checking BioBuddy enabled status:', error);
      this.bioBuddyEnabled = false;
    }

    await this.fetchComponentOptions();
    await this.generateWorkflow();
  },
  methods: {
    async fetchComponentOptions() {
      try {
        const workflow = await getWorkflow(this.workflowId);
        if (!workflow) {
          console.error("Failed to load workflow");
          return;
        }

        this.componentOptions = workflow.components.map(component => ({
          name: component.name,
          type: component.name, // Assuming type is the same as name, adjust if different
          inputs: component.input,
          outputs: component.output,
          description: this.getDescriptionString(component)
        }));
      } catch (error) {
        console.error("Error fetching component options:", error);
      }
    },
    getDescriptionString(component: ComponentModel_Output): string {
      // Assuming 'description' is a property within 'input', and converting it to string if it exists
      const inputKeys = Object.keys(component.input || {});
      const description = inputKeys.length > 0 ? component.input[inputKeys[0]].description : null;
      return description || "No description available";
    },
    async generateWorkflow() {
      try {
        const workflow = await getWorkflow(this.workflowId);
        if (!workflow) {
          console.error("Failed to load workflow");
          return;
        }

        const nodes: Node[] = [];
        const edges: Edge[] = [];

        // Create a map for component inputs and outputs
        const componentIOMap: { [key: string]: { inputs: string[], outputs: string[], type: string, description: string } } = {};
        workflow.components.forEach(component => {
          const inputs = Object.keys(component.input || {});
          const outputs = Object.keys(component.output || {});
          const description = this.getDescriptionString(component);
          componentIOMap[component.name] = { inputs, outputs, type: component.name, description }; // Assuming type is the same as name, adjust if different
        });

        // Process workflow_components to create nodes
        workflow.workflow_components.forEach(component => {
          const { inputs, outputs, type, description } = componentIOMap[component.name] || { inputs: [], outputs: [], type: '', description: '' };
          const nodeData: Node = {
            id: component.component_id,
            name: component.name,
            type: type, // Use the type from componentIOMap
            data: {
              description: description,
              inputs,
              outputs,
              jobIds: component.job_ids,
              draggable: false
            },
            position: { x: 100, y: 100 } // Default position, should be calculated based on layout logic
          };
          nodes.push(nodeData);
        });

        // Process mappings to create edges
        workflow.workflow_components.forEach(component => {
          component.mappings?.forEach(mapping => {
            edges.push({
              id: `e${mapping.source_component_id}-to-${component.component_id}`,
              source: mapping.source_component_id,
              target: component.component_id,
              sourceHandle: `${mapping.source_component_id}-output-${mapping.source_path[0]}`,
              targetHandle: `${component.component_id}-input-${mapping.target_path[0]}`,
            });
          });
        });

        this.elements = { nodes, edges };

        // Send the workflow update
        this.sendWorkflowUpdate();
      } catch (error) {
        console.error("Error generating workflow:", error);
      }
    },
    async sendWorkflowUpdate() {
      const workflowUpdate: WorkflowSchemaModel_Input = {
        workflow_id: this.workflowId,
        components: this.componentOptions.map(option => ({
          name: option.name,
          input: Object.keys(option.inputs).reduce((acc, key) => {
            acc[key] = option.inputs[key];
            return acc;
          }, {} as Record<string, PropertyModel_Output>),
          output: Object.keys(option.outputs).reduce((acc, key) => {
            acc[key] = option.outputs[key];
            return acc;
          }, {} as Record<string, PropertyModel_Output>)
        }) as ComponentModel_Input),
        workflow_components: this.elements.nodes.map(node => ({
          name: node.name,
          component_id: node.id,
          job_ids: node.data.jobIds,
          mappings: this.elements.edges
            .filter(edge => edge.target === node.id)
            .map(edge => ({
              source_path: [edge.sourceHandle?.split('-output-')[1]],
              target_path: [edge.targetHandle?.split('-input-')[1]],
              source_component_id: edge.sourceHandle?.split('-output-')[0],
              error: null
            }) as MappingModel),
          error: null,
          defaults: [],
          jobs_errors: []
        }) as WorkflowComponentModel),
        error: null,
        valid: true
      };
      try {
        await sendWorkflowUpdate(workflowUpdate);
        console.log('Workflow updated successfully');
      } catch (error) {
        console.error('Error updating workflow:', error);
      }
    },
    addComponent(option: any) {
      // Generate a new UUID for the node ID
      const newNodeId = uuidv4();

      // Create the new node
      const newNode: Node = {
        id: newNodeId,
        name: option.name,
        type: option.type, // Use type from option
        data: {
          description: option.description,
          inputs: Object.keys(option.inputs || {}),
          outputs: Object.keys(option.outputs || {}),
          jobIds: [],
          draggable: false
        },
        position: { x: 100, y: 100 }, // Default position
      };

      // Add the new node to the elements
      this.elements?.nodes.push(newNode);

      // Send the workflow update
      this.sendWorkflowUpdate();
    },
    onConnect(params: { source: string, target: string, sourceHandle: string, targetHandle: string }) {
      const newEdge: Edge = {
        id: `e${params.sourceHandle}-to-${params.targetHandle}`,
        source: params.sourceHandle.split('-output-')[0],
        target: params.targetHandle.split('-input-')[0],
        sourceHandle: params.sourceHandle,
        targetHandle: params.targetHandle,
      };
      this.elements.edges.push(newEdge);

      // Send the workflow update
      this.sendWorkflowUpdate();
    },
    onDeleteNode(nodeId: string) {
      const index = this.elements.nodes.findIndex(node => node.id === nodeId);
      if (index !== -1) {
        this.elements.nodes.splice(index, 1);
        this.elements.edges = this.elements.edges.filter(edge => edge.source !== nodeId && edge.target !== nodeId);

        // Send the workflow update
        this.sendWorkflowUpdate();
      }
    },
    openSettings(nodeId: string) {
      const node = this.elements.nodes.find(node => node.id === nodeId);
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
    onNodeDragStopHandler(event: { node: Node }) {
      const id = event.node.id;
      const newPosition = event.node.position;

      // Find the node with the given id and update its position
      const nodeToUpdate = this.elements.nodes.find((node: Node) => node.id === id);
      if (nodeToUpdate) {
        nodeToUpdate.position = newPosition;

        // Send the workflow update
        this.sendWorkflowUpdate();
      }
    },
    async onExperimentNameChange(newExperimentName: string) {
      const store = useDrugDiscoveryStore();
      await store.changeExperimentName(this.experiment.experimentId as string, newExperimentName);
    },
    openNodeDialog(nodeId: string) {
      const node = this.elements.nodes.find(node => node.id === nodeId);
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
