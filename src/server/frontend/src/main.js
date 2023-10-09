import 'bootswatch/dist/lux/bootstrap.min.css'; // Added this :boom:

import './assets/main.css';
import store from './storage';

import { createApp } from 'vue';
import App from './App.vue';
import router from './router';

const app = createApp(App);

app.use(router);
app.use(store);

app.mount('#app');
