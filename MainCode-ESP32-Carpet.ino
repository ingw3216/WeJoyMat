#include <FastLED.h>
#include <WiFi.h>
#include <WebServer.h>
#include <SPIFFS.h>
#include <ArduinoJson.h>
#include <esp_now.h>
#include <esp_wifi.h>
#include <HardwareSerial.h>
#include <DFRobotDFPlayerMini.h>
// WiFi设置
const char* ssid = "ESP32-MagicCarpet";
const char* password = "12345678";
WebServer server(80);
// 添加振动护腕的MAC地址
uint8_t wristbandMacAddress[] = {0xEC, 0xDA, 0x3B, 0xAA, 0x9F, 0x74}; 
// ESP-NOW通信设置
typedef struct vibration_message {
  bool activate;
  int duration; // 毫秒
} vibration_message;
// 新版本中使用esp_now_send_status_t作为第二个参数
void OnDataSent(const uint8_t *mac_addr, esp_now_send_status_t status) {
  Serial.print("发送状态: ");
  Serial.println(status == ESP_NOW_SEND_SUCCESS ? "成功" : "失败");
}
// DFPlayer音频播放器设置
HardwareSerial myDFSerial(2);  // 使用Serial2
DFRobotDFPlayerMini myDFPlayer;
bool dfPlayerReady = false;
// 游戏模式设置
#define MODE_VISUAL_ONLY 0           // 只有视觉
#define MODE_VISUAL_AUDIO 1          // 视觉+听觉
#define MODE_VISUAL_HAPTIC 2         // 视觉+触觉
#define MODE_ALL_SENSORY 3           // 视觉+听觉+触觉
volatile int currentGameMode = MODE_ALL_SENSORY;  // 默认为全感官模式
volatile int currentBrightness = 64;
// 引脚定义
#define LED_PIN     15
#define LED_COUNT   640                       // LED引脚15，总数192
const int parentSensorPins[] = {17, 5, 18, 19, 21};  
const int childSensorPins[] = {12, 14, 27, 26, 25};
const int motorPin = 32;                      // 震动电机引脚32

// 记录游戏数据
const int MAX_PLAYS = 200;                     // 每个游戏块最大记录游玩次数
int playCount[5] = {0, 0, 0, 0, 0};           // 每块板子的游玩次数
unsigned long reactionTimes[5][MAX_PLAYS];    // 存储每块板子的反应时间
unsigned long pressStart[5] = {0, 0, 0, 0, 0};   // 板子按下的时间点
unsigned long lightStart[5] = {0, 0, 0, 0, 0};   // 板子亮起的时间点

// 新增：按压序列记录数据结构
const int MAX_SEQUENCE_RECORDS = 1000;        // 最大序列记录数
struct SequenceRecord {
  int playCount;      // 第几次按压
  int blockNum;       // 按压的块号(1-5)
  unsigned long reactTime; // 反应时间
};
SequenceRecord sequenceData[MAX_SEQUENCE_RECORDS];
int totalSequenceCount = 0;                   // 总按压序列计数

// LED状态数组
CRGB leds[LED_COUNT];
volatile bool parentLedState[5] = {false, false, false, false, false};  // 父母侧LED初始状态为亮
volatile bool childLedState[5] = {false, false, false, false, false}; // 儿童侧LED初始状态为灭
unsigned long ledOffTime[5] = {0, 0, 0, 0, 0};  // 父母侧LED关闭的时间
const unsigned long ledResetDelay = 2000;       // 重新点亮LED的延时设置(2秒)

// 中断相关
volatile bool parentSensorTriggered[5] = {false, false, false, false, false};
volatile bool childSensorTriggered[5] = {false, false, false, false, false};
volatile bool ledNeedReset[5] = {false, false, false, false, false};
unsigned long lastDebounceTime[10] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
const unsigned long debounceDelay = 100;       // 100毫秒的防抖时间

// 用于蜂鸣器和振动电机的控制变量
bool motorActive = false;
unsigned long motorEndTime = 0;

// 记录反应时间
unsigned long reactStartTime[5] = {0, 0, 0, 0, 0}; 

// 定时器设置
hw_timer_t *timerParent = NULL;
hw_timer_t *timerChild = NULL;
portMUX_TYPE timerMux = portMUX_INITIALIZER_UNLOCKED;

// 父母侧传感器中断处理函数
void IRAM_ATTR parentSensorISR0() { if (!parentLedState[0]) parentSensorTriggered[0] = true; }
void IRAM_ATTR parentSensorISR1() { if (!parentLedState[1]) parentSensorTriggered[1] = true; }
void IRAM_ATTR parentSensorISR2() { if (!parentLedState[2]) parentSensorTriggered[2] = true; }
void IRAM_ATTR parentSensorISR3() { if (!parentLedState[3]) parentSensorTriggered[3] = true; }
void IRAM_ATTR parentSensorISR4() { if (!parentLedState[4]) parentSensorTriggered[4] = true; }

// 儿童侧传感器中断处理函数
void IRAM_ATTR childSensorISR0() { if (childLedState[0]) childSensorTriggered[0] = true; }
void IRAM_ATTR childSensorISR1() { if (childLedState[1]) childSensorTriggered[1] = true; }
void IRAM_ATTR childSensorISR2() { if (childLedState[2]) childSensorTriggered[2] = true; }
void IRAM_ATTR childSensorISR3() { if (childLedState[3]) childSensorTriggered[3] = true; }
void IRAM_ATTR childSensorISR4() { if (childLedState[4]) childSensorTriggered[4] = true; }

// 父母侧定时器中断处理函数
void IRAM_ATTR onParentTimer() {
  portENTER_CRITICAL_ISR(&timerMux);
  unsigned long currentTime = millis();
  
  // 检查是否有父母侧LED需要重新点亮
  for (int i = 0; i < 5; i++) {
    if (!parentLedState[i] && currentTime - ledOffTime[i] >= ledResetDelay) {
      ledNeedReset[i] = true;  // 标记需要重置LED
    }
  }
  portEXIT_CRITICAL_ISR(&timerMux);
}

// 儿童侧定时器中断处理函数
void IRAM_ATTR onChildTimer() {
  // 儿童侧不需要定时重置LED，由父母侧触发
}

void setup() {
  Serial.begin(115200);
  Serial.println("\n\n开始初始化...");

  // 设置FastLED
  FastLED.addLeds<WS2812, LED_PIN, GRB>(leds, LED_COUNT);
  FastLED.setBrightness(currentBrightness);  // 设置亮度为16

  // 设置引脚模式
  ledcAttach(motorPin, 5000, 8);
  myDFSerial.begin(9600, SERIAL_8N1, 22,23);  // DFPlayer连接
  Serial.println("初始化 DFPlayer Mini...");
  if (!myDFPlayer.begin(myDFSerial)) {
    Serial.println("DFPlayer初始化失败");
    dfPlayerReady = false;
  } else {
    Serial.println("DFPlayer初始化成功");
    myDFPlayer.volume(30);  // 设置默认音量（范围：0~30）
    dfPlayerReady = true;
  }

  for (int i = 0; i < 5; i++) {
    pinMode(parentSensorPins[i], INPUT);
    pinMode(childSensorPins[i], INPUT);
    // 清空反应时间数组
    for (int j = 0; j < MAX_PLAYS; j++) {
      reactionTimes[i][j] = 0;
    }
  }

  for (int i = 0; i < MAX_SEQUENCE_RECORDS; i++) {
    sequenceData[i].playCount = 0;
    sequenceData[i].blockNum = 0;
    sequenceData[i].reactTime = 0;
  }

  // 添加父母侧外部中断
  attachInterrupt(digitalPinToInterrupt(parentSensorPins[0]), parentSensorISR0, RISING);
  attachInterrupt(digitalPinToInterrupt(parentSensorPins[1]), parentSensorISR1, RISING);
  attachInterrupt(digitalPinToInterrupt(parentSensorPins[2]), parentSensorISR2, RISING);
  attachInterrupt(digitalPinToInterrupt(parentSensorPins[3]), parentSensorISR3, RISING);
  attachInterrupt(digitalPinToInterrupt(parentSensorPins[4]), parentSensorISR4, RISING);

  // 添加儿童侧外部中断
  attachInterrupt(digitalPinToInterrupt(childSensorPins[0]), childSensorISR0, RISING);
  attachInterrupt(digitalPinToInterrupt(childSensorPins[1]), childSensorISR1, RISING);
  attachInterrupt(digitalPinToInterrupt(childSensorPins[2]), childSensorISR2, RISING);
  attachInterrupt(digitalPinToInterrupt(childSensorPins[3]), childSensorISR3, RISING);
  attachInterrupt(digitalPinToInterrupt(childSensorPins[4]), childSensorISR4, RISING);

  // 设置定时器 - 父母
  timerParent = timerBegin(1000000);                 // 80分频，1us = 80个周期
  timerAttachInterrupt(timerParent, &onParentTimer); //边缘触发
  timerAlarm(timerParent, 100000, true, 0);      // 100ms检查一次

  // 设置定时器 - 儿童
  timerChild = timerBegin(1000000);           // 定时器1，80分频，向上计数
  timerAttachInterrupt(timerChild, &onChildTimer);   // 边缘触发
  timerAlarm(timerChild, 100000, true, 0);      // 100ms检查一次

  // 初始化WiFi和Web服务器
  setupWiFi();
  esp_wifi_set_channel(1, WIFI_SECOND_CHAN_NONE);
  WiFi.mode(WIFI_AP_STA); // 设置为AP+STA模式以同时支持Web服务器和ESP-NOW

  if (esp_now_init() != ESP_OK) {
    Serial.println("ESP-NOW初始化失败");
    return;
  }
  Serial.println("ESP-NOW初始化成功");

  // 注册发送回调
  esp_now_register_send_cb(OnDataSent);
  
  // 注册振动护腕为对等设备
  esp_now_peer_info_t peerInfo = {};
  memcpy(peerInfo.peer_addr, wristbandMacAddress, 6);
  peerInfo.channel = 1;  
  peerInfo.encrypt = false;
  
  // 添加对等设备前先检查是否已存在
  if (esp_now_is_peer_exist(wristbandMacAddress)) {
    esp_now_del_peer(wristbandMacAddress);
  }

  // 添加对等设备
  esp_err_t addStatus = esp_now_add_peer(&peerInfo);
  if (addStatus != ESP_OK) {
    Serial.print("无法添加对等设备，错误代码: ");
    Serial.println(addStatus);
    return;
  }
  
  Serial.println("对等设备添加成功");
  Serial.print("目标设备MAC: ");
  for(int i = 0; i < 6; i++) {
    Serial.printf("%02X", wristbandMacAddress[i]);
    if(i < 5) Serial.print(":");
  }
  Serial.println();

  // 初始化LED状态：父母侧亮，儿童侧灭
  for (int i = 1; i <= 10; i++) {
    ledTurnOff(i);  // 儿童侧灭
  }

  Serial.println("初始化成功，游戏开始!");
}

void loop() {
  unsigned long currentTime = millis();
  server.handleClient();  // 处理Web服务器请求

  // 处理父母侧传感器触发
  for (int i = 0; i < 5; i++) {
    if (parentSensorTriggered[i]) {
      if (currentTime - lastDebounceTime[i] > debounceDelay) {
        lastDebounceTime[i] = currentTime;
        
        // 确认父母侧传感器状态
        if (digitalRead(parentSensorPins[i]) == HIGH && !parentLedState[i]) {
          // 父母按下游戏块i
          ledTurnOn(i + 1);  // 关闭父母侧LED
          ledTurnOn(10 - i);  // 点亮儿童侧对应的LED (11-i)
          if (currentGameMode == MODE_VISUAL_AUDIO || currentGameMode == MODE_ALL_SENSORY) {
          playSoundForBlock(i + 1);
        }
          // 记录反应时间开始
          reactStartTime[i] = currentTime;
        }
      }
      
      portENTER_CRITICAL(&timerMux);
      parentSensorTriggered[i] = false;  // 重置触发标志
      portEXIT_CRITICAL(&timerMux);
    }
  }

  // 处理儿童侧传感器触发
  for (int i = 0; i < 5; i++) {
    if (childSensorTriggered[i]) {
      if (currentTime - lastDebounceTime[i + 5] > debounceDelay) {
        lastDebounceTime[i + 5] = currentTime;
        
        // 确认儿童侧传感器状态
        if (digitalRead(childSensorPins[i]) == HIGH && childLedState[i]) {
          // 儿童按下游戏块
          int parentIndex = 4 - i;  // 对应的父母侧索引 (反向映射)
          ledTurnOff(parentIndex + 1);  
          ledTurnOff(10 - parentIndex);
          // 计算反应时间并存储
          unsigned long reactionTime = currentTime - reactStartTime[parentIndex];
          if (playCount[parentIndex] < MAX_PLAYS) {
            reactionTimes[parentIndex][playCount[parentIndex]] = reactionTime;
            playCount[parentIndex]++;
          }
          if (totalSequenceCount < MAX_SEQUENCE_RECORDS) {
            sequenceData[totalSequenceCount].playCount = totalSequenceCount + 1;
            sequenceData[totalSequenceCount].blockNum = parentIndex + 1;  // 1-5的块号
            sequenceData[totalSequenceCount].reactTime = reactionTime;
            totalSequenceCount++;
          }
          if (currentGameMode == MODE_VISUAL_AUDIO || currentGameMode == MODE_ALL_SENSORY) {
            playSoundForBlock(parentIndex + 1);
          }
          if (currentGameMode == MODE_VISUAL_HAPTIC || currentGameMode == MODE_ALL_SENSORY) {
            startMotor();
          }   
          // 设置父母侧LED重置时间
          ledOffTime[parentIndex] = currentTime;
        }
      }  
      portENTER_CRITICAL(&timerMux);
      childSensorTriggered[i] = false;  // 重置触发标志
      portEXIT_CRITICAL(&timerMux);
    }
  }
}

void ledTurnOn(int ScreenNum) {
  int start = (ScreenNum - 1) * 64;
  int end = start + 63;
  CRGB color;
  
  // 根据屏幕编号设置不同颜色
  if (ScreenNum <= 5) {
    color = CRGB::Blue;  // 父母侧蓝色
  } else {
    color = CRGB::Green; // 儿童侧绿色
  }

  for (int i = start; i <= end; i++) {
    leds[i] = color;
  }
  // 更新LED状态
  if (ScreenNum <= 5) {
    parentLedState[ScreenNum - 1] = true;
  } else {
    childLedState[ScreenNum - 6] = true;
  }

  FastLED.show();
}
void ledTurnOff(int ScreenNum) {
  int start = (ScreenNum - 1) * 64;
  int end = start + 63;
  
  for (int i = start; i <= end; i++) {
    leds[i] = CRGB::Black;
  }
  // 更新LED状态
  if (ScreenNum <= 5) {
    parentLedState[ScreenNum - 1] = false;
  } else {
    childLedState[ScreenNum - 6] = false;
  }
  FastLED.show();
}
void playSoundForBlock(int blockNum) {
  if (!dfPlayerReady) return;
  int fileNumber = 0;
  switch (blockNum) {
    case 1: fileNumber = 1; break; // 播放0001.mp3
    case 2: fileNumber = 2; break; // 播放0002.mp3
    case 3: fileNumber = 3; break; // 播放0003.mp3
    case 4: fileNumber = 4; break; // 播放0005.mp3
    case 5: fileNumber = 5; break; // 播放0006.mp3
    default: return;
  }
  myDFPlayer.play(fileNumber);  // 播放对应文件，会立即打断之前播放
}

void startMotor() {
  // 保留本地蜂鸣器控制，移除本地振动电机控制
  if (currentGameMode == MODE_VISUAL_HAPTIC || currentGameMode == MODE_ALL_SENSORY) {
    // 发送振动命令到护腕
    vibration_message msg;
    msg.activate = true;
    msg.duration = 300; // 500ms振动
  
    // 使用esp_now_send发送数据
    esp_err_t result = esp_now_send(wristbandMacAddress, (uint8_t *)&msg, sizeof(msg));
  
    if (result == ESP_OK) {
      Serial.println("振动命令发送成功");
    } else {
      Serial.println("振动命令发送失败");
    }
  }
}
// 初始化WiFi
void setupWiFi() {
  // 配置固定IP地址
  IPAddress local_IP(192, 168, 4, 1);     // 固定IP地址
  IPAddress gateway(192, 168, 4, 1);      // 网关地址
  IPAddress subnet(255, 255, 255, 0);     // 子网掩码
  
  // 配置AP的IP地址
  WiFi.softAPConfig(local_IP, gateway, subnet);
  
  // 启动AP
  WiFi.softAP(ssid, password);

  // 设置Web服务器路由
  server.on("/", HTTP_GET, []() {
    String html = "<html><head><meta charset='UTF-8'><title>亲子互动游戏数据</title>";
    html += "<meta name='viewport' content='width=device-width, initial-scale=1'>";
    html += "<style>";
    html += "body{font-family:Arial;margin:20px;}";
    html += "h1,h2{color:#333;}";
    html += "table{border-collapse:collapse;width:100%;margin-bottom:20px;}";
    html += "th,td{border:1px solid #ddd;padding:8px;text-align:center;}";
    html += "th{background-color:#f2f2f2;}";
    html += "tr:nth-child(even){background-color:#f9f9f9;}";
    html += "select, button{font-size:18px;padding:10px 20px;margin-top:10px;}";
    html += ".tables-container{display:flex;gap:20px;flex-wrap:wrap;}";
    html += ".table-box{flex:1;min-width:300px;}";
    html += "</style></head><body>";

    html += "<h1>亲子互动游戏数据</h1>";

    html += "<h2>选择游玩模式</h2>";
    html += "<form action='/set_mode' method='POST'>";
    html += "<select name='mode'>";
    html += "<option value='0' " + String(currentGameMode == MODE_VISUAL_ONLY ? "selected" : "") + ">只有视觉</option>";
    html += "<option value='1' " + String(currentGameMode == MODE_VISUAL_AUDIO ? "selected" : "") + ">视觉+听觉</option>";
    html += "<option value='2' " + String(currentGameMode == MODE_VISUAL_HAPTIC ? "selected" : "") + ">视觉+触觉</option>";
    html += "<option value='3' " + String(currentGameMode == MODE_ALL_SENSORY ? "selected" : "") + ">视听触觉全部</option>";
    html += "</select>";
    html += "<button type='submit'>提交</button>";
    html += "</form>";

    html += "<h2>选择LED亮度</h2>";
    html += "<form action='/set_brightness' method='POST'>";
    html += "<select name='brightness'>";
    html += "<option value='64' " + String(currentBrightness == 64 ? "selected" : "") + ">低亮度</option>";
    html += "<option value='128' " + String(currentBrightness == 128 ? "selected" : "") + ">中亮度</option>";
    html += "<option value='200' " + String(currentBrightness == 200 ? "selected" : "") + ">高亮度</option>";
    html += "<option value='255' " + String(currentBrightness == 255 ? "selected" : "") + ">全亮度</option>";
    html += "</select>";
    html += "<button type='submit'>提交</button>";
    html += "</form>";

    html += "<form action='/reset' method='get'><button type='submit'>重置所有数据</button></form>";

    html += "<h2>统计信息</h2>";
    html += "<p>总按压次数: " + String(totalSequenceCount) + "</p>";
    html += "<p>各块按压次数: ";
    for (int i = 0; i < 5; i++) {
      html += "Block" + String(i+1) + "(" + String(playCount[i]) + "次) ";
    }
    html += "</p>";

    // 表格部分横向布局
    html += "<div class='tables-container'>";

    // 左表格：反应时间数据
    html += "<div class='table-box'>";
    html += "<h2>GameBlock反应时间(ms)</h2>";
    html += "<table><tr><th>playCount</th><th>block1</th><th>block2</th><th>block3</th><th>block4</th><th>block5</th></tr>";
    int maxPlayCount = 0;
    for (int i = 0; i < 5; i++) {
      if (playCount[i] > maxPlayCount) {
        maxPlayCount = playCount[i];
      }
}
    for (int j = 0; j < maxPlayCount; j++) {
      html += "<tr><td>" + String(j + 1) + "</td>";
      for (int i = 0; i < 5; i++) {
        if (j < playCount[i]) {
          html += "<td>" + String(reactionTimes[i][j]) + "</td>";
        } else {
          html += "<td>-</td>";
        }
      }
      html += "</tr>";
    }
    html += "</table></div>";

    // 右表格：按压序列数据
    html += "<div class='table-box'>";
    html += "<h2>游玩顺序反应时间</h2>";
    html += "<table><tr><th>playCount</th><th>blockNum</th><th>reactTime(ms)</th></tr>";
    for (int i = 0; i < totalSequenceCount; i++) {
      html += "<tr>";
      html += "<td>" + String(sequenceData[i].playCount) + "</td>";
      html += "<td>" + String(sequenceData[i].blockNum) + "</td>";
      html += "<td>" + String(sequenceData[i].reactTime) + "</td>";
      html += "</tr>";
    }
    html += "</table></div>";
    
    html += "</div>";  // end of tables-container

    html += "</body></html>";
    server.send(200, "text/html", html);
  });

  server.on("/set_mode", HTTP_POST, []() {
    String mode = server.arg("mode");
    currentGameMode = mode.toInt();  // 更新游戏模式
    server.sendHeader("Location", "/");
    server.send(303);
  });

  server.on("/set_brightness", HTTP_POST, []() {
    String bval = server.arg("brightness");
    currentBrightness = bval.toInt();  // 更新亮度
    FastLED.setBrightness(currentBrightness);  // 实时生效
    server.sendHeader("Location", "/");
    server.send(303);
  });

  server.on("/reset", HTTP_GET, []() {
    for (int i = 0; i < 5; i++) {
      playCount[i] = 0;  // 重置每块板子的游玩次数
    }
    totalSequenceCount = 0;
    for (int i = 0; i < MAX_SEQUENCE_RECORDS; i++) {
      sequenceData[i].playCount = 0;
      sequenceData[i].blockNum = 0;
      sequenceData[i].reactTime = 0;
    }
    server.sendHeader("Location", "/");
    server.send(303);
  });

  server.begin();
  Serial.println("Web服务器已启动");
}
