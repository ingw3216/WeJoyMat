#include <esp_now.h>
#include <WiFi.h>
#include <esp_wifi.h>

// 振动电机引脚定义
#define MOTOR_PIN 10 // 根据实际连接修改

// 定义消息结构
typedef struct vibration_message {
  bool activate;
  int duration; // 毫秒
} vibration_message;

// 当接收到数据时的回调函数
void OnDataRecv(const esp_now_recv_info *info, const uint8_t *incomingData, int len) {
  char macStr[18];
  snprintf(macStr, sizeof(macStr), "%02x:%02x:%02x:%02x:%02x:%02x",
           info->src_addr[0], info->src_addr[1], info->src_addr[2], 
           info->src_addr[3], info->src_addr[4], info->src_addr[5]);
  Serial.print("收到来自: ");
  Serial.println(macStr);
  
   // 处理接收到的数据
  if (len == sizeof(vibration_message)) {
    vibration_message msg;
    memcpy(&msg, incomingData, sizeof(msg));
  
    Serial.print("收到消息 - activate: ");
    Serial.print(msg.activate);
    Serial.print(", duration: ");
    Serial.println(msg.duration);

    if (msg.activate) {
    // 激活振动电机
      Serial.print("激活振动，持续时间: ");
      Serial.print(msg.duration);
      Serial.println("ms");

      ledcWrite(MOTOR_PIN, 255);
      delay(msg.duration);
      ledcWrite(MOTOR_PIN, 0);;

      Serial.println("振动结束");
    }else {
    Serial.println("收到的数据长度不匹配！");
    // 打印原始数据用于调试
    Serial.print("原始数据: ");
    for(int i = 0; i < len; i++) {
      Serial.printf("%02X ", incomingData[i]);
    }
    Serial.println();
  }
  }
}

void setup() {
  Serial.begin(115200);
  ledcAttach(MOTOR_PIN, 5000, 8);
  Serial.println("开始初始化XIAO-ESP32-C3...");

  // 初始化WiFi为站点模式
  WiFi.mode(WIFI_STA);
  // 设置WiFi信道为1（与主ESP32保持一致）
  esp_wifi_set_channel(1, WIFI_SECOND_CHAN_NONE);

  // 打印本机MAC地址（用于主ESP32设置）
  Serial.print("护腕MAC地址: ");
  Serial.println(WiFi.macAddress());
 
  // 打印WiFi信道信息
  uint8_t channel;
  wifi_second_chan_t second;
  esp_wifi_get_channel(&channel, &second);
  Serial.print("WiFi信道: ");
  Serial.println(channel);

  // 初始化ESP-NOW
  if (esp_now_init() != ESP_OK) {
    Serial.println("ESP-NOW初始化失败");
    return;
  }
  Serial.println("ESP-NOW初始化成功");  

  // 注册接收回调
  //esp_now_register_recv_cb(OnDataRecv);
  esp_err_t result = esp_now_register_recv_cb(OnDataRecv);
  if (result != ESP_OK) {
    Serial.print("注册接收回调失败，错误代码: ");
    Serial.println(result);
    return;
  }
  Serial.println("接收回调注册成功");


  // 启动测试 - 短暂振动表示准备就绪
  ledcWrite(MOTOR_PIN, 255);
  delay(500);
  ledcWrite(MOTOR_PIN, 0);
  Serial.println("ESP32-C3 振动护腕初始化完成，等待指令...");
}

void loop() {
  // 添加心跳信息，每5秒打印一次状态
  static unsigned long lastHeartbeat = 0;
  unsigned long currentTime = millis();
  
  if (currentTime - lastHeartbeat > 5000) {
    lastHeartbeat = currentTime;
    Serial.println("设备运行中，等待ESP-NOW消息...");
  }
  
  delay(10);
}