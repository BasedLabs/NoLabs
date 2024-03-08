<template>
    <q-list>
      <q-item-label header>BioBuddy Chat</q-item-label>
      <q-separator />
      <q-item v-for="(message, index) in messages" :key="index">
        <q-item-section>
          <div v-if="message.type === 'function'" class="function-message">
            <q-item-label class="text-h7 q-mb-sm">
              <div class="q-pb-sm q-pt-sm text-bold">{{ displayName(message.role) }}</div>
              <p> <q-icon name="check_circle" color="purple"></q-icon>
                <span class="text-h7 text-purple q-ml-sm"> {{ displayContent(message) }} </span>
              </p>
              <ul>
                Params:
                <li v-for="(param, index) in message.message.parameters" :key="index">
                  {{ param.name }}: {{ param.value }}
                </li>
              </ul>
            </q-item-label>
          </div>
          <div v-else>
            <q-item-label class="text-h7 q-mb-sm">
              <div class="q-pb-sm q-pt-sm text-bold">{{ displayName(message.role) }}</div>
              <p>{{ displayContent(message) }}</p>
            </q-item-label>
          </div>
        </q-item-section>
      </q-item>
    </q-list>
    <div class="q-pa-md">
      <q-input v-model="newMessage" label="Type a message..." dense filled>
        <template v-slot:append>
          <q-btn icon="send" flat @click="sendMessage" :disable="sending">
            <q-spinner size="20px" v-if="sending"></q-spinner>
          </q-btn>
        </template>
      </q-input>
    </div>
</template>

<script lang="ts">
import { QList, QItem, QItemLabel, QSeparator, QItemSection, QInput, QBtn, QSpinner } from 'quasar';
import { defineComponent } from 'vue';
import { loadConversationApi, sendMessageApi } from 'src/features/biobuddy/api';
import {FunctionCall, Message, type RegularMessage} from "src/api/client";

export interface FunctionMapping {
  name: string;
  function: (parameters: any) => void;
}

export default defineComponent({
  name: 'BioBuddyChat',
  components: {
    QList, QItem, QItemLabel, QSeparator, QItemSection, QInput, QBtn, QSpinner
  },
  props: {
    experimentId: {
      type: String,
      required: true,
    },
    functionMappings: {
      type: Array as () => FunctionMapping[],
      required: false,
    }
  },
  data() {
    return {
      drawer: false,
      messages: [] as Message[],
      newMessage: '',
      sending: false,
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
      const response = await sendMessageApi(this.experimentId, this.newMessage);
      const newMessageResponse = response.biobuddy_response as Message;
      this.messages.push(newMessageResponse);

      if (newMessageResponse.type === 'function') {
        const functionCall = newMessageResponse.message as FunctionCall;
        this.invokeFunction(functionCall.function_name);
      }
      this.newMessage = '';
      await this.loadConversation();
      this.sending = false;
    },
    invokeFunction(functionName: string) {
      const mapping = this.functionMappings?.find(m => m.name === functionName);
      if (mapping && typeof mapping.function === 'function') {
        // Accepts storage functions
        mapping.function(this.experimentId);
      }
    },
    displayContent(message: Message) {
      if (message.type === 'text') {
        const response_message = message.message as RegularMessage;
        return response_message.content;
      }
      const response_message = message.message as FunctionCall;
      return response_message.function_name;
    },
    displayName(role: string) {
      return role === 'user' ? 'You' : 'Biobuddy';
    }
  },
  mounted() {
    this.loadConversation();
  },
});
</script>
