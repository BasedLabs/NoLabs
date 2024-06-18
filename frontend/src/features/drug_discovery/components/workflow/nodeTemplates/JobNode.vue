<template>
  <q-card dark bordered v-if="!nodeData">
    <q-spinner color="primary" size="3em" />
  </q-card>
  <q-card dark bordered v-if="nodeData" style="max-width: 400px">
    <q-card-section v-if="nodeData?.data.error">
      <div class="row no-wrap items-center">
        <q-space />
        <q-btn v-if="nodeData?.data.error" color="red" class="q-pm-md">
          <q-icon left name="warning" />
          <div>{{ nodeData?.data.error }}</div>
        </q-btn>
      </div>
    </q-card-section>
    <q-card-section>
      <JobNodeContent :nodeId="nodeId" :name="nodeData?.name" :description="nodeData?.description"/>
    </q-card-section>
    <q-card-section>
      <NodeHandles :nodeId="nodeId" :inputs="nodeData?.data.inputs" :outputs="nodeData?.data.outputs" />
    </q-card-section>
    <q-card-actions align="around">
      <q-btn v-if="!isRunning" @click="openDialogue" color="positive" icon="play_arrow" label="Start" />
      <q-btn v-else>
        <q-spinner color="primary" size="20px" />
        Running...
      </q-btn>
      <q-btn @click="openDialogue" label="Extended view" />
      <q-btn @click="deleteNode" color="negative" icon="delete" label="Delete" />
    </q-card-actions>
  </q-card>
</template>

<script>
import JobNodeContent from './JobNodeContent.vue';
import NodeHandles from './NodeHandles.vue';
import { useWorkflowStore } from 'src/features/drug_discovery/components/workflow/storage';

export default {
  name: 'JobNode',
  components: {
    JobNodeContent,
    NodeHandles
  },
  props: {
    nodeId: {
      type: String,
      required: true
    },
    onDeleteNode: Function,
    onOpenSettings: Function,
    onOpenDialog: Function
  },
  data() {
    return {
      nodeData: null
    };
  },
  computed: {
    isRunning() {
      const workflowStore = useWorkflowStore();
      return workflowStore.runningComponentIds.includes(this.nodeId);
    }
  },
  watch: {
    'useWorkflowStore().runningComponentIds': function(newVal, oldVal) {
      if (newVal !== oldVal) {
        this.updateRunningStatus();
      }
    }
  },
  mounted() {
    const workflowStore = useWorkflowStore();
    this.nodeData = workflowStore.getNodeById(this.nodeId);
    console.log(this.nodeData.name);
  },
  methods: {
    deleteNode() {
      if (this.onDeleteNode) {
        this.onDeleteNode(this.nodeId);
      }
    },
    openSettings() {
      if (this.onOpenSettings) {
        this.onOpenSettings(this.nodeId);
      }
    },
    openDialogue() {
      if (this.onOpenDialog) {
        this.onOpenDialog(this.nodeId);
      }
    },
    updateRunningStatus() {
      const workflowStore = useWorkflowStore();
      this.isRunning = workflowStore.runningComponentIds.includes(this.nodeId);
    }
  }
};
</script>
