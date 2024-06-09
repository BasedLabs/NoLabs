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
                {{ displayName(message.role) }}</div>
              <div v-for="(functionCall, fcIndex) in message.message" :key="`function-${fcIndex}`">
              <p> <q-icon name="check_circle" color="purple"></q-icon>
                <span class="text-h7 text-purple q-ml-sm"> {{ functionCall.function_name }} </span>
              </p>
              <ul>
                Params:
                <li v-for="(param, pIndex) in functionCall.arguments" :key="`param-${fcIndex}-${pIndex}`">
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
                <div class="markdown-content" v-html="renderMarkdown(message.message.content)"></div>
                <q-btn flat v-if="message.role === 'user'" size="sm" icon="edit" @click="startEditing(index)"></q-btn>
                <q-btn flat v-if="isLastUserMessageWithoutResponse(index) && !sending" label="Regenerate" icon-right="autorenew" @click="sendQuery(message.message.content)"></q-btn>
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
import {Notify} from 'quasar';
import { defineComponent } from 'vue';
import {editMessageApi, loadConversationApi, saveMessageApi, sendQueryApi} from 'src/features/biobuddy/api';
import {FunctionCall, type FunctionParam, Message} from "src/refinedApi/client";
import MarkdownIt from 'markdown-it';
import {useBioBuddyStore} from "./storage";
import {obtainErrorResponse} from "../../api/errorWrapper";

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
    renderMarkdown(text: string) {
      const md = new MarkdownIt();
      let result = md.render(text);
      result = result.replace(
        /<h1>/g,
        '<h1 style="font-size: 1.5em; margin-top: 0.5em; margin-bottom: 0.5em; line-height: 1.2; font-weight: bold;">'
      );
      result = result.replace(/<h2>/g, '<h2 style="font-size: 1.3em; margin-top: 0.3em; margin-bottom: 0.3em;">');
      result = result.replace(/<h3>/g, '<h3 style="font-size: 1.1em; margin-top: 0.3em; margin-bottom: 0.3em;">');
      return result;
    },
    async sendMessage() {
      if (!this.experimentId || !this.newMessage.trim()) return;

      try {
        const response = await saveMessageApi(this.experimentId, this.newMessage);
        const savedMessage = response.saved_message as Message;

        this.messages.push(savedMessage);
        const queryContent = this.newMessage;
        this.newMessage = '';
        await this.sendQuery(queryContent);
      } catch (error) {
        console.error("Failed to send or process message:", error);
      }
    },
    isLastUserMessageWithoutResponse(index: number) {
      return index === this.messages.length - 1 && this.messages[index].role === 'user';
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
        await this.sendQuery(this.editMessage);
      } catch (error) {
        console.error("Failed to edit message:", error);
      }
      this.editIndex = -1;
      this.editMessage = '';
    },
    async sendQuery(message: string) {
      this.sending = true;
      const queryResponse = await sendQueryApi(this.experimentId, message);
      const errorResponse = obtainErrorResponse(queryResponse);
      if (errorResponse) {
        this.sending = false;
        for (const error of errorResponse.errors) {
          Notify.create({
            type: "negative",
            closeBtn: 'Close',
            message: error
          });
        }
        return;
      }
      const newMessageResponse = queryResponse.biobuddy_response as Message;
      this.messages.push(newMessageResponse);
      if (newMessageResponse.type === 'function') {
        const functionCall = newMessageResponse.message as FunctionCall[];
        await this.invokeFunctions(functionCall);
      }

      this.sending = false;
    },
    async invokeFunctions(functionCalls: FunctionCall[]) {
      for (const functionCall of functionCalls) {
        const mapping = this.functionMappings?.find(m => m.name === functionCall.function_name);
        if (mapping && typeof mapping.function === 'function') {
          await mapping.function(functionCall.data);
        }
      }
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
        return (message.message as FunctionCall[]).map(fc => fc.arguments ? fc.arguments : []);
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
</style>
