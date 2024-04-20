<template>
  <q-drawer show-if-above
                :width="drawerWidth"
                    >
    <q-list>
      <q-item-label header>BioBuddy Chat</q-item-label>
      <q-separator />
      <q-item v-for="(message, index) in messages" :key="index">
        <q-item-section>
          <div v-if="message.type === 'function'" class="function-message">
            <q-item-label class="text-h7 q-mb-sm">
              <div class="q-pb-sm q-pt-sm text-bold">
                <img src="/Biobuddy_icon.svg" class="custom-icon" />
                {{ displayName(message.role) }}</div>
              <div v-for="(functionCall, fcIndex) in message.message" :key="`function-${fcIndex}`">
              <p> <q-icon name="check_circle" color="purple"></q-icon>
                <span class="text-h7 text-purple q-ml-sm"> {{ displayContent(message) }} </span>
              </p>
              <ul>
                Params:
                <li v-for="(param, pIndex) in functionCall.parameters" :key="`param-${fcIndex}-${pIndex}`">
                  {{ param.name }}: {{ param.value }}
                </li>
              </ul>
              </div>
            </q-item-label>
            </div>
          <div v-else>
            <q-item-label class="text-h7 q-mb-sm">
              <div v-if="editIndex !== index">
                <div class="q-pb-sm q-pt-sm text-bold">{{ displayName(message.role) }}</div>
                <p>{{ displayContent(message) }}</p>
                <q-btn flat v-if="message.role === 'user'" label="Edit" @click="startEditing(index)"></q-btn>
              </div>
              <div v-else>
                <q-input v-model="editMessage" dense filled></q-input>
                <q-btn flat label="Save and Submit" @click="saveEdit(index)"></q-btn>
                <q-btn flat label="Cancel" @click="cancelEdit"></q-btn>
              </div>
            </q-item-label>
          </div>
        </q-item-section>
      </q-item>
    </q-list>
    <div class="q-pa-md">
      <q-input v-model="newMessage"  label="Type a message..." dense filled @keyup.enter="sendMessage">
        <template v-slot:append>
          <q-btn icon="send" flat @click="sendMessage" :disable="sending">
            <q-spinner size="20px" v-if="sending"></q-spinner>
          </q-btn>
        </template>
      </q-input>
    </div>


    <div
      class="drawer-resize-handle"
      @mousedown="startResizing"
    ></div>
  </q-drawer>
</template>

<script lang="ts">
import {QList, QItem, QItemLabel, QSeparator, QItemSection, QInput, QBtn, QSpinner} from 'quasar';
import { defineComponent } from 'vue';
import {editMessageApi, loadConversationApi, saveMessageApi, sendQueryApi} from 'src/features/biobuddy/api';
import {FunctionCall, type FunctionParam, Message, type RegularMessage} from "src/api/client";
import {useBioBuddyStore} from "./storage";

export interface FunctionMapping {
  name: string;
  function: (parameters: any) => void;
}

export default defineComponent({
  name: 'BioBuddyChat',
  props: {
    experimentId: {
      type: String,
      required: true,
    }
  },
  data() {
    return {
      drawer: true,
      messages: [] as Message[],
      newMessage: '',
      sending: false,
      functionMappings: [] as FunctionMapping[],
      drawerWidth: 500, // Initial width of the drawer
      isResizing: false,
      initialMouseX: 0,
      editIndex: -1,
      editMessage: '',
    };
  },
  methods: {
    async loadConversation() {
      if (!this.experimentId) return;
      const response = await loadConversationApi(this.experimentId);
      this.messages = response.messages;
    },
    async sendMessage() {
      if (!this.experimentId || !this.newMessage.trim()) return;
      this.sending = true;

      try {
        const response = await saveMessageApi(this.experimentId, this.newMessage);
        const savedMessage = response.saved_message as Message;

        this.messages.push(savedMessage);
        const queryContent = this.newMessage;
        this.newMessage = '';
        const queryResponse = await sendQueryApi(this.experimentId, queryContent);
        const newMessageResponse = queryResponse.biobuddy_response as Message;
        this.messages.push(newMessageResponse);
        if (newMessageResponse.type === 'function') {
          const functionCall = newMessageResponse.message as FunctionCall[];
          await this.invokeFunctions(functionCall);
        }

      } catch (error) {
        console.error("Failed to send or process message:", error);
      }
      this.sending = false;
    },
    startEditing(index: number) {
      this.editIndex = index;
      this.editMessage = this.messages[index].message.content;
    },
    cancelEdit() {
      this.editIndex = -1;
      this.editMessage = '';
    },
    async saveEdit(index: number) {
      if (!this.experimentId || !this.editMessage.trim()) return;
      const messageId = this.messages[index].id; // Assuming each message has a unique ID
      try {
        await editMessageApi(this.experimentId, messageId, this.editMessage);
        this.messages[index].message.content = this.editMessage;
        this.messages = this.messages.slice(0, index + 1); // Remove messages after the edited one

        this.sending = true;
        const queryResponse = await sendQueryApi(this.experimentId, this.editMessage);
        const newMessageResponse = queryResponse.biobuddy_response as Message;
        this.messages.push(newMessageResponse);
        if (newMessageResponse.type === 'function') {
          const functionCall = newMessageResponse.message as FunctionCall[];
          await this.invokeFunctions(functionCall);
        }

        this.sending = false;

      } catch (error) {
        console.error("Failed to edit message:", error);
      }
      this.editIndex = -1;
      this.editMessage = '';
    },
    async invokeFunctions(functionCalls: FunctionCall[]) {
      for (const functionCall of functionCalls) {
        const mapping = this.functionMappings?.find(m => m.name === functionCall.function_name);
        if (mapping && typeof mapping.function === 'function') {
          await mapping.function(functionCall.data);
        }
      }
    },
    displayContent(message: Message) {
      if (message.type === 'text') {
        const response_message = message.message as RegularMessage;
        return response_message.content;
      }
      const functionCalls = message.message as FunctionCall[];
      return functionCalls.map(fc => fc.function_name).join(', ');
    },
    displayName(role: string) {
      return role === 'user' ? 'You' : 'Biobuddy';
    },
    startResizing(event: MouseEvent) {
      this.isResizing = true;
      this.initialMouseX = event.clientX;
      window.addEventListener('mousemove', this.resizeDrawer);
      window.addEventListener('mouseup', this.stopResizing);
    },
    resizeDrawer(event: MouseEvent) {
      if (!this.isResizing) return;
      const deltaX = event.clientX - this.initialMouseX;
      this.drawerWidth += deltaX;
      this.initialMouseX = event.clientX; // Update initial X to the current position
    },
    stopResizing() {
      this.isResizing = false;
      window.removeEventListener('mousemove', this.resizeDrawer);
      window.removeEventListener('mouseup', this.stopResizing);
    },
    getParameters(message: Message): Array<Array<FunctionParam>> {
      if (message.type === 'function') {
        return (message.message as FunctionCall[]).map(fc => fc.parameters ? fc.parameters : []);
      }
      return [];
    }
  },
  async mounted() {
    await this.loadConversation();
    const bioBuddyStore = useBioBuddyStore();
    this.functionMappings = [
      { name: 'query_rcsb_pdb_by_id', function: bioBuddyStore.invokeQueryRcsbPdbEventHandlers },
      { name: 'query_rcsb_pdb_by_protein_names', function: bioBuddyStore.invokeQueryRcsbPdbEventHandlers },
      { name: 'query_chembl', function: bioBuddyStore.invokeQueryChemblEventHandlers },
      { name: 'query_chembl_by_condition', function: bioBuddyStore.invokeQueryChemblEventHandlers }
    ];
  },
});
</script>


<style scoped>
.drawer-resize-handle {
  cursor: ew-resize;
  position: absolute;
  top: 0;
  right: 0; /* Adjust based on your layout, ensuring it's reachable for resizing */
  width: 20px;
  height: 100%;
}

.custom-icon {
  width: 15px; /* Example size */
  height: 15px; /* Example size */
  vertical-align: top; /* Aligns icon with text if necessary */
}
</style>
